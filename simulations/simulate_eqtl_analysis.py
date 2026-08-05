import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink
import gzip
import subprocess
import tempfile





def load_in_cis_variants(gene_cis_variant_file):
	f = open(gene_cis_variant_file)
	ordered_variant_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ordered_variant_names.append(data[1])
	f.close()

	return np.asarray(ordered_variant_names)


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

def standardize_geno(geno_mat):
	std_geno_mat = np.copy(geno_mat)
	n_snps = std_geno_mat.shape[0]

	for snp_iter in range(n_snps):
		snp_mean = np.mean(std_geno_mat[snp_iter,:])
		snp_std = np.std(std_geno_mat[snp_iter,:], ddof=1)
		if snp_std == 0.0 or np.isfinite(snp_std) == False:
			std_geno_mat[snp_iter,:] = 0.0
		else:
			std_geno_mat[snp_iter,:] = (std_geno_mat[snp_iter,:] - snp_mean)/snp_std

	return std_geno_mat


def create_mapping_from_gene_id_to_causal_effects(causal_variant_gene_effect_size_file):
	f = gzip.open(causal_variant_gene_effect_size_file,'rt')
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
		effect = float(data[6])
		if gene_id not in mapping:
			mapping[gene_id] = []
		mapping[gene_id].append((var_id, effect))
	f.close()
	return mapping

def marginal_snp_effects_and_se(gene_expression, geno_mat):
	"""
	Compute marginal SNP effect sizes and standard errors from univariate
	linear regressions with an intercept.

	Parameters
	----------
	gene_expression : ndarray, shape (n_samples,)
		Expression vector.
	geno_mat : ndarray, shape (n_snps, n_samples)
		Genotype matrix.

	Returns
	-------
	betas : ndarray, shape (n_snps,)
		Marginal OLS effect size for each SNP.
	ses : ndarray, shape (n_snps,)
		Standard error of each marginal effect size.
	"""
	y = np.asarray(gene_expression, dtype=float)
	G = np.asarray(geno_mat, dtype=float)

	if y.ndim != 1:
		raise ValueError("gene_expression must be 1D")
	if G.ndim != 2:
		raise ValueError("geno_mat must be 2D")
	if G.shape[1] != y.shape[0]:
		raise ValueError("geno_mat must have shape (n_snps, n_samples)")

	n = y.shape[0]
	if n <= 2:
		raise ValueError("Need at least 3 samples")

	# Center y and each SNP genotype vector to account for intercept
	y_c = y - np.mean(y)
	G_c = G - np.mean(G, axis=1, keepdims=True)

	# Sxx_j = sum_i (g_ji - mean_g_j)^2
	sxx = np.sum(G_c**2, axis=1)

	# Identify monomorphic / zero-variance SNPs
	valid = sxx > 0

	betas = np.full(G.shape[0], np.nan)
	ses = np.full(G.shape[0], np.nan)

	if not np.any(valid):
		return betas, ses

	# beta_j = sum_i Gc_ji * yc_i / sum_i Gc_ji^2
	numer = G_c[valid] @ y_c
	betas_valid = numer / sxx[valid]

	# Residual sum of squares for each SNP:
	# RSS_j = sum(y_c^2) - beta_j^2 * Sxx_j
	# because with intercept and centered variables, fitted values are beta_j * G_cj
	y_ss = np.sum(y_c**2)
	rss = y_ss - (betas_valid**2) * sxx[valid]

	# Numerical guard against tiny negative values from floating point
	rss = np.maximum(rss, 0.0)

	sigma2 = rss / (n - 2)
	ses_valid = np.sqrt(sigma2 / sxx[valid])

	betas[valid] = betas_valid
	ses[valid] = ses_valid

	return betas, ses


def run_susie_fine_mapping(gene_expression, std_geno_mat, cis_variant_names, gene_name, susie_r_script, susie_num_effects=10):
	"""
	Run SuSiE fine-mapping with individual-level expression and standardized genotypes.

	std_geno_mat has shape (n_snps, n_samples). susieR expects X with shape
	(n_samples, n_snps), so the matrix is transposed before writing.
	"""
	valid_snps = np.all(np.isfinite(std_geno_mat), axis=1)
	if np.sum(valid_snps) == 0:
		return ''

	with tempfile.TemporaryDirectory() as tmp_dir:
		x_file = os.path.join(tmp_dir, 'susie_X.txt')
		y_file = os.path.join(tmp_dir, 'susie_y.txt')
		variant_file = os.path.join(tmp_dir, 'susie_variants.txt')
		output_file = os.path.join(tmp_dir, 'susie_results.txt')

		np.savetxt(x_file, np.transpose(std_geno_mat[valid_snps, :]), delimiter='\t')
		np.savetxt(y_file, gene_expression, delimiter='\t')
		np.savetxt(variant_file, cis_variant_names[valid_snps], fmt='%s', delimiter='\t')

		cmd = [
			'Rscript',
			susie_r_script,
			x_file,
			y_file,
			variant_file,
			gene_name,
			output_file,
			str(susie_num_effects)
		]
		try:
			subprocess.run(cmd, check=True, capture_output=True, text=True)
		except subprocess.CalledProcessError as error:
			print('SuSiE fine-mapping failed for ' + gene_name)
			print(error.stderr)
			raise

		with open(output_file) as f:
			return f.read()


def create_mapping_from_variant_id_to_alleles(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
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

######################
# Command line args
######################
causal_variant_gene_effect_size_file = sys.argv[1]
est_eqtl_effect_size_file = sys.argv[2]
gene_ld_summary_file = sys.argv[3]
onek_genomes_plink_filestem = sys.argv[4]
eqtl_sample_size = int(sys.argv[5])
simulation_iter = int(sys.argv[6])
individual_expression_file = sys.argv[7]
susie_fine_mapping_file = sys.argv[8] if len(sys.argv) > 8 else None
genotype_sample_mapping_file = sys.argv[9] if len(sys.argv) > 9 else None
susie_r_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'run_susie_fine_mapping.R')



# set random seed
np.random.seed(simulation_iter)

# Create mapping from gene id to vector of causal effects
gene_id_to_causal_effects = create_mapping_from_gene_id_to_causal_effects(causal_variant_gene_effect_size_file)



# Open output file handle
t = gzip.open(est_eqtl_effect_size_file, 'wt')
t.write('gene\tvariant\tchr\tsnp_pos\tA0\tA1\teqtl_effect_size\teqtl_effect_size_se\tN\tmaf\tgenotype_sdev\n')

# Open output file handle
t_indi = gzip.open(individual_expression_file, 'wt')
t_indi.write('gene\tindividual\texpression\n')

# Open SuSiE fine-mapping output file handle
if susie_fine_mapping_file is not None:
	t_susie = gzip.open(susie_fine_mapping_file, 'wt')
	t_susie.write('gene\tvariant\tsusie_pip\tsusie_posterior_mean\tsusie_posterior_sd\tsusie_credible_sets\tsusie_cs_min_abs_corr\n')
else:
	t_susie = None

# Loop through chromosomes
for chrom_num in range(1,23):
	print(chrom_num)
	# string of chromosome name
	chrom_string = 'chr' + str(chrom_num)

	# Load in chromosome plink data
	(bim, fam, G) = read_plink(onek_genomes_plink_filestem + str(chrom_num))

	if eqtl_sample_size > G.shape[1]:
		raise ValueError('Requested eqtl_sample_size=' + str(eqtl_sample_size) + ' but only ' + str(G.shape[1]) + ' genotype samples are available.')

	if chrom_num == 1:
		# Randomly sample individuals for eqtl analysis
		sample_indices = np.sort(np.random.choice(G.shape[1], size=eqtl_sample_size, replace=False))
		if genotype_sample_mapping_file is not None:
			np.savetxt(genotype_sample_mapping_file, sample_indices, fmt='%d')

	# Create mapping from variant id to index
	rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))

	# Create mapping from rsid to a0, a1
	rsid_to_alleles = create_mapping_from_variant_id_to_alleles(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))


	# Loop through genes / filter to only genes on this chromosome
	f = open(gene_ld_summary_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_chrom_num = data[1]
		if line_chrom_num != chrom_string:
			continue
		gene_name = data[0]
		counter = counter + 1
		if gene_name not in gene_id_to_causal_effects:
			continue
		gene_cis_variant_file = data[4]
		n_cis_variants = int(data[3])
		if n_cis_variants < 10:
			continue
		cis_variant_names = load_in_cis_variants(gene_cis_variant_file)
		print(gene_name)
		# Extract variant indices corresponding to cis_variant names
		cis_variant_indices = []
		cis_variant_a0s = []
		cis_variant_a1s = []
		cis_variant_chroms = []
		cis_variant_poss = []
		for variant_name in cis_variant_names:
			cis_variant_indices.append(rsid_to_genotype_index[variant_name])
			var_infos = rsid_to_alleles[variant_name]
			cis_variant_a0s.append(var_infos[0])
			cis_variant_a1s.append(var_infos[1])
			cis_variant_chroms.append(var_infos[2])
			cis_variant_poss.append(var_infos[3])
		cis_variant_indices = np.asarray(cis_variant_indices)
		cis_variant_a0s = np.asarray(cis_variant_a0s)
		cis_variant_a1s = np.asarray(cis_variant_a1s)
		cis_variant_chroms = np.asarray(cis_variant_chroms)
		cis_variant_poss = np.asarray(cis_variant_poss)

		# Extract genotype matrix
		geno_mat = G[cis_variant_indices,:].compute()
		geno_mat = geno_mat[:, sample_indices]
		n_samp = geno_mat.shape[1]

		# Allele frequency within the eqtl sample, so it matches the genotype scale the
		# marginal effect sizes below are estimated on
		allele_freqs = np.mean(geno_mat, axis=1)/2.0
		mafs = np.minimum(allele_freqs, 1.0 - allele_freqs)

		# Standardize genotype matrix
		std_geno_mat = standardize_geno(geno_mat)

		# Standardized genotype matrix
		genotype_sdevs = np.std(geno_mat, axis=1, ddof=1)

		# extract causal effects for gene
		causal_effect_tupler = gene_id_to_causal_effects[gene_name]
		tmp_variant_names = []
		tmp_causal_effects = []
		for tuply in causal_effect_tupler:
			tmp_variant_names.append(tuply[0])
			tmp_causal_effects.append(tuply[1])
		tmp_variant_names = np.asarray(tmp_variant_names)
		tmp_causal_effects = np.asarray(tmp_causal_effects)

		# QUick error check
		if np.array_equal(cis_variant_names, tmp_variant_names) == False:
			print('assumption erororo')
			pdb.set_trace()

		# Get genetic component of gene expression
		genetic_gene = np.dot(np.transpose(std_geno_mat), tmp_causal_effects)

		# Simulate gene expression
		gene_var = np.var(genetic_gene,ddof=1)
		if gene_var > 1:
			print(gene_var)
			resid_var = 1e-5
		else:
			resid_var = 1.0 - gene_var
		gene_expression = np.random.normal(loc=genetic_gene, scale=np.sqrt(resid_var))

		# Run eqtl mapping
		eqtl_beta, eqtl_beta_se = marginal_snp_effects_and_se(gene_expression, geno_mat)

		# Run SuSiE fine-mapping
		if t_susie is not None:
			susie_results = run_susie_fine_mapping(gene_expression, std_geno_mat, cis_variant_names, gene_name, susie_r_script)
			t_susie.write(susie_results)
			t_susie.flush()


		# Print to output file
		for var_iter, variant_id in enumerate(cis_variant_names):
			t.write(gene_name + '\t' + variant_id + '\t' + str(cis_variant_chroms[var_iter]) + '\t' + str(cis_variant_poss[var_iter]) + '\t' + cis_variant_a0s[var_iter] + '\t' + cis_variant_a1s[var_iter] + '\t' + str(eqtl_beta[var_iter]) + '\t' + str(eqtl_beta_se[var_iter]) + '\t' + str(n_samp) + '\t' + str(mafs[var_iter]) + '\t' + str(genotype_sdevs[var_iter]) + '\n')

		# Print to individual output file
		for indi_iter, indi_expr in enumerate(gene_expression):
			t_indi.write(gene_name + '\t' + str(indi_iter) + '\t' + str(indi_expr) + '\n')


	f.close()
t.close()
t_indi.close()
if t_susie is not None:
	t_susie.close()
print(est_eqtl_effect_size_file)
if susie_fine_mapping_file is not None:
	print(susie_fine_mapping_file)
