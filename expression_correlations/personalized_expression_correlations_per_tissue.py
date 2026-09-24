import numpy as np
import sys
import pdb
import gzip
import argparse
from pandas_plink import read_plink
from scipy import stats
from scipy.optimize import minimize_scalar


# Bin edges of the 'borzoi_finer_magnitude_bins' annotation (must match annotate_variant_gene_pairs.py).
# Bins are assigned on |unstandardized borzoi effect| and are right-open: bins[k] <= value < bins[k+1].
finer_borzoi_magnitude_bins = [0.0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, np.inf]
finer_borzoi_magnitude_anno_name = 'borzoi_finer_magnitude_bins'
finer_borzoi_magnitude_category_prefix = 'magnitude_bin'

# Which point estimate to pull from the S-LDMC bootstrap_stats file
sldmc_estimate_column = 'mean'

# Number of Monte Carlo draws of the true causal effects used to compute the directional FSR
n_fsr_samples = 1000

# Minimum number of usable cis variants for a gene to be analyzed
min_cis_variants = 10


def extract_bin_index(value, bins):
	for bin_iter in range(len(bins)-1):
		if value >= bins[bin_iter] and value < bins[bin_iter+1]:
			return bin_iter
	return -1


def load_in_sldmc_finer_magnitude_bin_estimates(sldmc_summary_file, anno_name, category_prefix, n_bins, estimate_column):
	# Pull the per-magnitude-bin calibration slope, per-snp eQTL h2 and correlation from the
	# S-LDMC (cross-tissue meta-analyzed) bootstrap_stats file.
	# Columns (tab-separated, with header): annotation_name  category_name  output_name  mean  bootstrapped_mean  ...
	slopes = np.full(n_bins, np.nan)
	per_snp_h2s = np.full(n_bins, np.nan)
	corrs = np.full(n_bins, np.nan)

	f = open(sldmc_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			if estimate_column not in data:
				print('assumption eroror: ' + estimate_column + ' column not found in ' + sldmc_summary_file)
				pdb.set_trace()
			estimate_index = data.index(estimate_column)
			continue
		if data[0] != anno_name:
			continue
		category_name = data[1]
		output_name = data[2]
		if category_name.startswith(category_prefix) == False:
			print('assumption eroror: unexpected category ' + category_name)
			pdb.set_trace()
		bin_index = int(category_name.split(category_prefix)[1])
		if bin_index < 0 or bin_index >= n_bins:
			print('assumption eroror: unexpected bin index ' + category_name)
			pdb.set_trace()
		estimate = float(data[estimate_index])
		if output_name == 'calibration_slope':
			slopes[bin_index] = estimate
		elif output_name == 'per_snp_eqtl_h2':
			per_snp_h2s[bin_index] = estimate
		elif output_name == 'correlation':
			corrs[bin_index] = estimate
	f.close()

	if np.any(np.isnan(slopes)) or np.any(np.isnan(per_snp_h2s)):
		print('assumption eroror: missing S-LDMC estimates for annotation ' + anno_name + ' in ' + sldmc_summary_file)
		pdb.set_trace()

	# Residual variance of the true causal effect given the (rescaled) borzoi prediction, per bin:
	# Var(beta | borzoi) = h2_b * (1 - corr_b^2). Fall back to h2_b if the correlation is undefined,
	# and clip at zero so a noisy bin never yields a negative variance.
	resid_vars = np.zeros(n_bins)
	for bin_iter in range(n_bins):
		if np.isfinite(corrs[bin_iter]):
			resid_vars[bin_iter] = per_snp_h2s[bin_iter]*(1.0 - np.square(corrs[bin_iter]))
		else:
			resid_vars[bin_iter] = per_snp_h2s[bin_iter]
	resid_vars[resid_vars < 0.0] = 0.0

	return slopes, per_snp_h2s, corrs, resid_vars


def create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file):
	# Columns (tab-separated, with header): gene  variant  chr  snp_pos  a0  a1  borzoi_effect_size
	f = gzip.open(est_borzoi_effect_size_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		var_id = data[1]
		chrom_num = data[2]
		snp_pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		effect = float(data[6])

		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()

		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, effect)
	f.close()
	return mapping


def create_mapping_from_gene_id_to_expression_vector(expr_file):
	# Residualized expression bed file: 4 leading columns (chr, start, end, gene_id) then one column per sample
	dicti = {}
	head_count = 0
	f = open(expr_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			expr_sample_names = np.asarray(data[4:])
			continue
		gene_id = data[3].split('.')[0]
		expr = np.asarray(data[4:]).astype(float)
		if gene_id in dicti:
			print('assumption erororr')
			pdb.set_trace()
		dicti[gene_id] = expr
	f.close()
	return dicti, expr_sample_names


def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
	mapping = {}
	n_snps = len(ordered_snps)
	for snp_iter in range(n_snps):
		snp_name = ordered_snps[snp_iter]
		if snp_name in mapping:
			print('asssumption erororo')
			pdb.set_trace()
		mapping[snp_name] = snp_iter
	return mapping


def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
	if len(snp_array) != len(a0_arr):
		print('assumption eorroro')
		pdb.set_trace()
	if len(snp_array) != len(a1_arr):
		print('assumption eorroro')
		pdb.set_trace()

	dicti = {}
	for ii, snp_id in enumerate(snp_array):
		if snp_id in dicti:
			print('assumpationoenroer')
			pdb.set_trace()
		dicti[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
	return dicti


def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num


def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_borzoi_effects):
	# Variants with a borzoi effect that are in the genotype data with matching alleles
	unique_vars = np.unique([*var_to_est_borzoi_effects])
	final_vars = []
	for var in unique_vars:
		if var not in rsid_to_genotype_index:
			continue
		geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])
		borzoi_alleles = var_to_est_borzoi_effects[var][4:6]
		if set(geno_alleles) != set(borzoi_alleles):
			continue
		final_vars.append(var)
	return np.asarray(final_vars)


def load_in_snp_gene_data(ordered_cis_variants, var_to_est_effects):
	effects = []
	alleles = []
	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_effects:
			print('assumption erororr')
			pdb.set_trace()
		var_info = var_to_est_effects[variant_id]
		effects.append(var_info[6])
		alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles)


def assign_finer_magnitude_bins(borzoi_effects_unstandardized, bins):
	bin_indices = []
	for effect in borzoi_effects_unstandardized:
		bin_index = extract_bin_index(np.abs(effect), bins)
		if bin_index == -1:
			print('assumption eroror: borzoi effect falls outside magnitude bins')
			pdb.set_trace()
		bin_indices.append(bin_index)
	return np.asarray(bin_indices)


def estimate_cis_snp_heritability_with_lrt(genotype_mat, expr_vec):
	# Cis-SNP heritability via maximum likelihood on a single variance component model
	# (y ~ N(0, h2*GRM + (1-h2)*I) up to a scale), with a likelihood ratio test against h2 = 0.
	# The p-value uses the 50:50 mixture of chi2(0) and chi2(1) since h2 = 0 is on the boundary.
	X = np.asarray(genotype_mat, dtype=float)
	y = np.asarray(expr_vec, dtype=float)
	y = y - np.mean(y)

	X = X - np.mean(X, axis=0)
	snp_sdevs = np.std(X, axis=0)
	valid_snps = np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
	X = X[:, valid_snps]
	snp_sdevs = snp_sdevs[valid_snps]
	if X.shape[1] == 0:
		return np.nan, np.nan

	X = X/snp_sdevs[None, :]
	n_samples = X.shape[0]
	grm = np.dot(X, np.transpose(X))/X.shape[1]
	eigenvalues, eigenvectors = np.linalg.eigh(grm)
	transformed_y = np.dot(np.transpose(eigenvectors), y)

	def log_likelihood(h2):
		variance_scale = h2*eigenvalues + (1.0 - h2)
		if np.any(variance_scale <= 0.0):
			return -np.inf
		residual_var = np.mean(np.square(transformed_y)/variance_scale)
		if residual_var <= 0.0:
			return -np.inf
		return -0.5*(n_samples*np.log(2.0*np.pi) + n_samples*np.log(residual_var) + np.sum(np.log(variance_scale)) + n_samples)

	null_log_likelihood = log_likelihood(0.0)
	opt = minimize_scalar(lambda h2: -log_likelihood(h2), bounds=(0.0, 0.999999), method='bounded')
	h2 = opt.x
	alt_log_likelihood = -opt.fun
	lrt_stat = np.maximum(2.0*(alt_log_likelihood - null_log_likelihood), 0.0)
	lrt_pvalue = 0.5*stats.chi2.sf(lrt_stat, df=1)
	return h2, lrt_pvalue


def compute_cis_grm(genotype_mat):
	# Standardize genotypes in sample (discarding variants with no variance) and return the cis GRM X X^T / m
	X = np.asarray(genotype_mat, dtype=float)
	X = X - np.mean(X, axis=0)
	snp_sdevs = np.std(X, axis=0)
	valid_snps = np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
	if np.sum(valid_snps) == 0:
		return None
	X = X[:, valid_snps]/snp_sdevs[valid_snps][None, :]
	return np.dot(X, np.transpose(X))/X.shape[1]


def estimate_cis_snp_heritability_with_he_regression(genotype_mat, expr_vec):
	# Haseman-Elston regression estimate of cis-SNP heritability with an analytic standard error.
	# Regress y_i y_k on GRM_ik (with an intercept) over all pairs i != k of the mean-centered expression y; the slope
	# divided by Var(y) is h2. Unlike the bounded MLE above this is unbiased under the single-component model and can be
	# negative, which matters when averaging over low-h2 genes (e.g. when comparing against predicted_cis_snp_h2).
	# Standard error: the estimate is a ratio of quadratic forms in y, R = (y^T W y) / (y^T y / n), with W symmetric and
	# zero-diagonal. Var(R) follows from Var(y^T M y) = 2 tr(M S M S), Cov(y^T M y, y^T P y) = 2 tr(M S P S) and the delta
	# method, with S the covariance of the centered y under the fitted single-component model h2*GRM + (1-h2)*I (h2
	# clipped to [0, 1]). Because W has zero diagonal this needs no kurtosis assumption on the residuals, only that they
	# are independent and homoscedastic. Resampling individuals is deliberately not used: under h2 = 0, y^T W y is a
	# degenerate U-statistic, for which the delete-one jackknife and the naive bootstrap over individuals overstate the
	# variance (~2x; checked by simulation).
	y = np.asarray(expr_vec, dtype=float)
	y = y - np.mean(y)
	n = len(y)
	grm = compute_cis_grm(genotype_mat)
	var_y = np.mean(np.square(y))
	if grm is None or n < 3 or var_y <= 0.0:
		return np.nan, np.nan

	# W_ik = (GRM_ik - mean off-diagonal GRM) / sum_{i != k} (GRM_ik - mean)^2 for i != k, and 0 on the diagonal,
	# so that y^T W y is the OLS slope (with intercept) of y_i y_k on GRM_ik over pairs i != k
	off_diag = ~np.eye(n, dtype=bool)
	W = grm - np.mean(grm[off_diag])
	W[~off_diag] = 0.0
	sxx = np.sum(np.square(W))
	if sxx <= 0.0:
		return np.nan, np.nan
	W = W/sxx
	he_h2 = float(np.dot(y, np.dot(W, y))/var_y)

	# Analytic standard error. Covariance of the centered y under the fitted model:
	# C (s2_g GRM + s2_e I) C = s2_g GRM + s2_e C, since the in-sample GRM is already double-centered (GRM 1 = 0)
	h2_fit = np.clip(he_h2, 0.0, 1.0)
	centering = np.eye(n) - 1.0/n
	Sigma = var_y*(h2_fit*grm + (1.0 - h2_fit)*centering)
	W_Sigma = np.dot(W, Sigma)
	var_q1 = 2.0*np.sum(W_Sigma*np.transpose(W_Sigma))   # Var(y^T W y)            = 2 tr(W S W S)
	cov_q1_q2 = 2.0*np.sum(W_Sigma*Sigma)/n             # Cov(y^T W y, y^T y / n) = 2 tr(W S S) / n
	var_q2 = 2.0*np.sum(np.square(Sigma))/np.square(n)  # Var(y^T y / n)          = 2 tr(S S) / n^2
	var_h2 = (var_q1 - 2.0*he_h2*cov_q1_q2 + np.square(he_h2)*var_q2)/np.square(var_y)
	he_h2_se = float(np.sqrt(np.maximum(var_h2, 0.0)))
	return he_h2, he_h2_se


def compute_directional_fsr(genotype_mat, rescaled_pred_expr, per_snp_mu, per_snp_sd, n_samp):
	# Monte Carlo directional false sign rate:
	# probability that the reported score (X mu) points the wrong way relative to a draw of the
	# true genetic expression (X beta^(s)), with beta_j^(s) ~ N(mu_j, sd_j^2) independently per snp.
	sampled_betas = np.random.normal(loc=per_snp_mu[None, :], scale=per_snp_sd[None, :], size=(n_samp, len(per_snp_mu)))
	# n_individuals X n_samples
	sampled_pred_expr = np.dot(genotype_mat, np.transpose(sampled_betas))

	# Center over individuals so the alignment corresponds to a covariance
	reported_score_centered = rescaled_pred_expr - np.mean(rescaled_pred_expr)
	sampled_pred_expr_centered = sampled_pred_expr - np.mean(sampled_pred_expr, axis=0)

	# Alignment between reported score and each sampled true score: (X mu)^T (X beta^(s))
	alignment_samples = np.dot(reported_score_centered, sampled_pred_expr_centered)

	directional_fsr = np.mean(alignment_samples < 0)
	return directional_fsr


def compute_predicted_cis_snp_heritability(rescaled_pred_expr, per_snp_sd):
	# Predicted (expected) cis-SNP heritability of expression given the rescaled borzoi predictions:
	# E[Var_i(X beta) | borzoi] with beta_j ~ N(mu_j, sd_j^2) independently per snp and X standardized in sample
	#   = Var_i(X mu) + sum_j sd_j^2
	# This is the closed form of the mean over Monte Carlo draws of Var_i(X beta^(s)); the LD cross-terms are
	# carried entirely by Var_i(X mu) because the residuals are independent across snps and diag(X^T X / n) = 1.
	# On the standardized-expression scale of the S-LDMC estimates, so directly comparable to cis_snp_h2.
	return np.var(rescaled_pred_expr) + np.sum(np.square(per_snp_sd))


def compute_per_bin_mean_genotype_variance(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, plink_genotype_stem, bins, chunk_size=10000):
	# Mean in-sample genotype variance (2p(1-p) under HWE; the square of the sdev used to standardize)
	# across variant-gene pairs falling in each finer magnitude bin, computed from this tissue's genotype
	# data. Uses the same genes, variants and filters as the main loop (genes with expression,
	# allele-matched variants, no missing calls, nonzero variance).
	n_bins = len(bins) - 1
	bin_sums = np.zeros(n_bins)
	bin_counts = np.zeros(n_bins)

	for chrom_num in range(1,23):
		print(chrom_num)
		(bim, fam, G) = read_plink(plink_genotype_stem + str(chrom_num))
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))

		# For each variant on this chromosome, the bins of its variant-gene pairs
		variant_to_bin_indices = {}
		for gene_id in [*gene_id_to_est_borzoi_effects]:
			gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
			if str(gene_chrom_num) != str(chrom_num):
				continue
			if gene_id not in gene_id_to_expression_vector:
				continue
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])
			if len(ordered_cis_variants) < min_cis_variants:
				continue
			for cis_variant in ordered_cis_variants:
				bin_index = extract_bin_index(np.abs(gene_id_to_est_borzoi_effects[gene_id][cis_variant][6]), bins)
				if bin_index == -1:
					print('assumption eroror: borzoi effect falls outside magnitude bins')
					pdb.set_trace()
				if cis_variant not in variant_to_bin_indices:
					variant_to_bin_indices[cis_variant] = []
				variant_to_bin_indices[cis_variant].append(bin_index)

		unique_variants = np.asarray([*variant_to_bin_indices])
		if len(unique_variants) == 0:
			continue
		unique_genotype_indices = np.asarray([rsid_to_genotype_index[variant] for variant in unique_variants])
		# Sort by genotype index (contiguous plink reads)
		ordering = np.argsort(unique_genotype_indices)
		unique_variants = unique_variants[ordering]
		unique_genotype_indices = unique_genotype_indices[ordering]

		# Genotype variance of each variant, in chunks to bound memory
		for chunk_start in range(0, len(unique_variants), chunk_size):
			chunk_variants = unique_variants[chunk_start:(chunk_start + chunk_size)]
			chunk_indices = unique_genotype_indices[chunk_start:(chunk_start + chunk_size)]
			geno_mat = (G[chunk_indices,:].compute())[:, genotype_sample_indices]
			observed_variants = np.any(np.isnan(geno_mat), axis=1) == False
			geno_vars = np.nanvar(geno_mat, axis=1)
			valid_variants = observed_variants & np.isfinite(geno_vars) & (geno_vars > 0.0)
			for ii, variant in enumerate(chunk_variants):
				if valid_variants[ii] == False:
					continue
				for bin_index in variant_to_bin_indices[variant]:
					bin_sums[bin_index] = bin_sums[bin_index] + geno_vars[ii]
					bin_counts[bin_index] = bin_counts[bin_index] + 1

	bin_mean_geno_vars = np.full(n_bins, np.nan)
	for bin_iter in range(n_bins):
		if bin_counts[bin_iter] > 0:
			bin_mean_geno_vars[bin_iter] = bin_sums[bin_iter]/bin_counts[bin_iter]
	return bin_mean_geno_vars, bin_counts


def run_expression_correlations(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, plink_genotype_stem, bin_slopes, bin_resid_vars, bin_tau2s, output_file):
	# Initialize output file
	t = open(output_file,'w')
	t.write('gene_id\traw_expression_correlation\trescaled_expression_correlation\texpression_FSR\texpression_FSR_af_specific\tcis_snp_h2\tcis_snp_h2_pvalue\tcis_snp_h2_he\tcis_snp_h2_he_se\tpredicted_cis_snp_h2\tpredicted_cis_snp_h2_af_specific\n')

	n_genes_analyzed = 0

	# Loop through chromsomes
	for chrom_num in range(1,23):
		print(chrom_num)

		##################################
		# Load in per-chrom-genotype data
		##################################
		# Load in chromosome plink data
		(bim, fam, G) = read_plink(plink_genotype_stem + str(chrom_num))
		# Create mapping from variant id to index
		rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
		# Create mapping from rsid to a0, a1
		rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))


		##################################
		# Loop through genes on this chromosome
		# (Analysis done seperately for each gene)
		##################################
		for gene_id in [*gene_id_to_est_borzoi_effects]:

			# Limit to genes on this chromosome
			gene_chrom_num = extract_gene_chrom_num(gene_id_to_est_borzoi_effects[gene_id])
			if str(gene_chrom_num) != str(chrom_num):
				continue

			# Gene needs both borzoi effects AND expression
			if gene_id not in gene_id_to_expression_vector:
				continue

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id])
			if len(ordered_cis_variants) < min_cis_variants:
				continue

			# Load in borzoi effects for gene (unstandardized: per allele)
			borzoi_effects_unstandardized, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])

			# Genotype indices of cis variants
			cis_genotype_indices = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				geno_alleles = rsid_to_snp_info[cis_variant][:2]
				# Genotype dosage below is (2 - G), which counts copies of bim a0.
				# Borzoi effects are per copy of their a1 allele, so flip when borzoi a1 is not bim a0.
				if borzoi_variant_alleles[var_index,:][1] != geno_alleles[0]:
					borzoi_effects_unstandardized[var_index] = -1.0*borzoi_effects_unstandardized[var_index]
			cis_genotype_indices = np.asarray(cis_genotype_indices)

			# Extract genotype matrix (variants X expression samples)
			geno_mat = (G[cis_genotype_indices,:].compute())[:, genotype_sample_indices]

			# drop_missing: discard any variant with a missing genotype call (matches S-LDMC default)
			observed_variants = np.any(np.isnan(geno_mat), axis=1) == False

			geno_mat = 2.0 - geno_mat

			# Standardize genotype (in sample); discard variants with no variance
			snp_means = np.nanmean(geno_mat, axis=1)
			snp_sdevs = np.nanstd(geno_mat, axis=1)
			valid_snps = observed_variants & np.isfinite(snp_sdevs) & (snp_sdevs > 0.0)
			if np.sum(valid_snps) < min_cis_variants:
				continue

			genotype_mat = np.transpose((geno_mat[valid_snps, :] - snp_means[valid_snps, None])/snp_sdevs[valid_snps, None])
			borzoi_vec_unstandardized = borzoi_effects_unstandardized[valid_snps]
			# Borzoi effects on the standardized genotype scale (the scale S-LDMC slopes were estimated on)
			borzoi_vec = borzoi_vec_unstandardized*snp_sdevs[valid_snps]

			expr_vec = gene_id_to_expression_vector[gene_id]
			if len(expr_vec) != genotype_mat.shape[0]:
				print('assumption eroror: expression and genotype sample mismatch')
				pdb.set_trace()

			# Raw prediction: X beta_borzoi
			raw_pred_expr = np.dot(genotype_mat, borzoi_vec)
			raw_corry = np.corrcoef(expr_vec, raw_pred_expr)[0,1]

			# Rescaled prediction: X (slope_bin * beta_borzoi), bins assigned on |unstandardized borzoi|
			bin_indices = assign_finer_magnitude_bins(borzoi_vec_unstandardized, finer_borzoi_magnitude_bins)
			per_snp_mu = bin_slopes[bin_indices]*borzoi_vec
			rescaled_pred_expr = np.dot(genotype_mat, per_snp_mu)
			rescaled_corry = np.corrcoef(expr_vec, rescaled_pred_expr)[0,1]

			# Directional FSR of the rescaled prediction
			# (a) residual variance constant within bin on the standardized scale: Var(eps_j) = h2_c(1 - corr_c^2)
			per_snp_sd = np.sqrt(bin_resid_vars[bin_indices])
			expression_fsr = compute_directional_fsr(genotype_mat, rescaled_pred_expr, per_snp_mu, per_snp_sd, n_fsr_samples)

			# (b) allele-frequency-specific residual variance (Supplementary Note eq. 35): Var(eps_j) = 2p_j(1-p_j) tau^2_c
			per_snp_sd_af_specific = snp_sdevs[valid_snps]*np.sqrt(bin_tau2s[bin_indices])
			if np.all(np.isfinite(per_snp_sd_af_specific)):
				expression_fsr_af_specific = compute_directional_fsr(genotype_mat, rescaled_pred_expr, per_snp_mu, per_snp_sd_af_specific, n_fsr_samples)
			else:
				expression_fsr_af_specific = np.nan

			# Predicted cis-SNP heritability given the rescaled borzoi predictions, under each residual-variance model
			predicted_cis_snp_h2 = compute_predicted_cis_snp_heritability(rescaled_pred_expr, per_snp_sd)
			if np.all(np.isfinite(per_snp_sd_af_specific)):
				predicted_cis_snp_h2_af_specific = compute_predicted_cis_snp_heritability(rescaled_pred_expr, per_snp_sd_af_specific)
			else:
				predicted_cis_snp_h2_af_specific = np.nan

			# Cis-SNP heritability of expression (and LRT p-value)
			cis_snp_h2, cis_snp_h2_pvalue = estimate_cis_snp_heritability_with_lrt(genotype_mat, expr_vec)
			# Haseman-Elston regression estimate (unbounded, can be negative) with analytic standard error
			cis_snp_h2_he, cis_snp_h2_he_se = estimate_cis_snp_heritability_with_he_regression(genotype_mat, expr_vec)

			t.write(gene_id + '\t' + str(raw_corry) + '\t' + str(rescaled_corry) + '\t' + str(expression_fsr) + '\t' + str(expression_fsr_af_specific) + '\t' + str(cis_snp_h2) + '\t' + str(cis_snp_h2_pvalue) + '\t' + str(cis_snp_h2_he) + '\t' + str(cis_snp_h2_he_se) + '\t' + str(predicted_cis_snp_h2) + '\t' + str(predicted_cis_snp_h2_af_specific) + '\n')
			t.flush()
			n_genes_analyzed = n_genes_analyzed + 1

	t.close()
	print(str(n_genes_analyzed) + ' genes analyzed')
	return





#####################
# Command line args
#####################
parser = argparse.ArgumentParser(description='Compute correlations between personalized (Borzoi-predicted) expression and observed expression in a single tissue.')
parser.add_argument('--sldmc-summary-file', dest='sldmc_summary_file', required=True, help='SLDMC cross-tissue meta-analyzed summary statistics file.')
parser.add_argument('--borzoi-results-file', dest='borzoi_results_file', required=True, help='Borzoi predicted variant effect sizes file for this tissue/sample (gzipped).')
parser.add_argument('--expression-file', dest='expression_file', required=True, help='Residualized, renormalized expression bed file for this tissue.')
parser.add_argument('--plink-genotype-stem', dest='plink_genotype_stem', required=True, help='Stem of chromosome-specific plink genotype files (chromosome number and .bed/.bim/.fam are appended).')
parser.add_argument('--genotype-sample-mapping-file', dest='genotype_sample_mapping_file', required=True, help='File mapping genotype sample indices to expression samples for this tissue.')
parser.add_argument('--expression-correlation-output-file', dest='expression_correlation_output_file', required=True, help='Output file containing personalized expression prediction correlations.')
args = parser.parse_args()

sldmc_summary_file = args.sldmc_summary_file
borzoi_results_file = args.borzoi_results_file
expression_file = args.expression_file
plink_genotype_stem = args.plink_genotype_stem
genotype_sample_mapping_file = args.genotype_sample_mapping_file
expression_correlation_output_file = args.expression_correlation_output_file


np.random.seed(1)

###########################
# Load in data
###########################

# Per finer-magnitude-bin S-LDMC estimates (calibration slope used to rescale borzoi effects;
# residual variance used to sample true causal effects for the FSR)
n_bins = len(finer_borzoi_magnitude_bins) - 1
bin_slopes, bin_per_snp_h2s, bin_corrs, bin_resid_vars = load_in_sldmc_finer_magnitude_bin_estimates(sldmc_summary_file, finer_borzoi_magnitude_anno_name, finer_borzoi_magnitude_category_prefix, n_bins, sldmc_estimate_column)
for bin_iter in range(n_bins):
	print(finer_borzoi_magnitude_category_prefix + str(bin_iter) + ' [' + str(finer_borzoi_magnitude_bins[bin_iter]) + ', ' + str(finer_borzoi_magnitude_bins[bin_iter+1]) + '): slope=' + str(bin_slopes[bin_iter]) + ' per_snp_h2=' + str(bin_per_snp_h2s[bin_iter]) + ' corr=' + str(bin_corrs[bin_iter]) + ' resid_var=' + str(bin_resid_vars[bin_iter]))


# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(borzoi_results_file)

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

# Create mapping from gene id to expression vector
gene_id_to_expression_vector, expr_sample_names = create_mapping_from_gene_id_to_expression_vector(expression_file)
if len(expr_sample_names) != len(genotype_sample_indices):
	print('assumption eroror: expression samples and genotype sample mapping have different lengths')
	pdb.set_trace()


# Allelic residual variance per bin (Supplementary Note eq. 34):
# tau^2_c = (Var_c(beta^s) - mu_c^2 Var_c(delta^s)) / mean_c(2p(1-p))
# The numerator is bin_resid_vars; the denominator is the mean in-sample genotype variance over the
# variant-gene pairs in the bin, computed from this tissue's genotype data.
print('Computing per-bin mean genotype variance')
bin_mean_geno_vars, bin_pair_counts = compute_per_bin_mean_genotype_variance(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, plink_genotype_stem, finer_borzoi_magnitude_bins)
bin_tau2s = bin_resid_vars/bin_mean_geno_vars
for bin_iter in range(n_bins):
	print(finer_borzoi_magnitude_category_prefix + str(bin_iter) + ': n_pairs=' + str(int(bin_pair_counts[bin_iter])) + ' mean_geno_var=' + str(bin_mean_geno_vars[bin_iter]) + ' tau2=' + str(bin_tau2s[bin_iter]))


###########################
# Run expression correlation analysis
###########################
run_expression_correlations(gene_id_to_est_borzoi_effects, genotype_sample_indices, gene_id_to_expression_vector, plink_genotype_stem, bin_slopes, bin_resid_vars, bin_tau2s, expression_correlation_output_file)
