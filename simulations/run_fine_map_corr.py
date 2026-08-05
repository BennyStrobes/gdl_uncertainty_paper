import numpy as np
import os
import sys
import pdb
import gzip
import pickle
from pandas_plink import read_plink
import time





def create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file):
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


def create_mapping_from_gene_id_to_variant_gene_annotations(sim_variant_gene_annotation_file):
	f = gzip.open(sim_variant_gene_annotation_file,'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[6:])
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
		anno = np.asarray(data[6:]).astype(float)
		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, anno)
	f.close()
	return mapping, anno_names


def create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file):	
	f = gzip.open(est_eqtl_effect_size_file,'rt')
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
		pos = data[3]
		a0 = data[4]
		a1 = data[5]
		if a0 == a1:
			print('assumption eroroor')
			pdb.set_trace()
		effect = float(data[6])
		se = float(data[7])
		eqtl_sample_size = float(data[8])
		desired_se = np.sqrt(1.0/eqtl_sample_size)
		if gene_id not in mapping:
			mapping[gene_id] = {}
		zed = effect/se
		std_effect = zed*desired_se
		if var_id in mapping[gene_id]:
			print('repeat snp error')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, pos, a0, a1, std_effect, desired_se)
	f.close()
	return mapping

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

def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_borzoi_effects, var_to_est_eqtl_effects):
	unique_vars = np.unique(np.hstack(([*var_to_est_borzoi_effects],[*var_to_est_eqtl_effects])))
	final_vars = []
	for var in unique_vars:
		if var not in rsid_to_genotype_index:
			continue
		geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])

		passing = True
		if var in var_to_est_borzoi_effects:
			borzoi_alleles = var_to_est_borzoi_effects[var][4:6]
			if set(geno_alleles) != set(borzoi_alleles):
				passing = False
		if var in var_to_est_eqtl_effects:
			eqtl_alleles = var_to_est_eqtl_effects[var][4:6]
			if set(geno_alleles) != set(eqtl_alleles):
				passing = False

		if passing == False:
			continue
		final_vars.append(var)
	return np.asarray(final_vars)


def load_in_snp_gene_eqtl_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []
	effect_ses = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.nan)
			effect_ses.append(np.nan)
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			effect_ses.append(var_info[7])
			alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles), np.asarray(effect_ses)

def load_in_snp_gene_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			effects.append(np.nan)
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_est_eqtl_effects[variant_id]
			effects.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.asarray(effects), np.asarray(alleles)

def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num

def generate_gene_ld_means(gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, gene_id_to_variant_gene_anno, onek_genomes_plink_filestem, anno_names):
	# Initialize output object
	gene_to_ld_means = {}
	
	# Loop through chromsomes
	for chrom_num in range(1,23):
		print(chrom_num)

		##################################
		# Load in per-chrom-genotype data
		##################################
		# string of chromosome name
		chrom_string = 'chr' + str(chrom_num)
		# Load in chromosome plink data
		(bim, fam, G) = read_plink(onek_genomes_plink_filestem + str(chrom_num))
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

			# Gene needs both borzoi effects AND eQTLs
			if gene_id not in gene_id_to_est_eqtl_effects:
				continue
			if gene_id not in gene_id_to_variant_gene_anno:
				continue

			# Extract ordered list of variants
			ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, gene_id_to_est_borzoi_effects[gene_id], gene_id_to_est_eqtl_effects[gene_id])
			# Sip genes with fewer than 10 variants
			if len(ordered_cis_variants) < 10:
				continue

			# Load in data for gene
			# eQTL
			eqtl_effects, eqtl_variant_alleles, eqtl_effect_ses = load_in_snp_gene_eqtl_data(ordered_cis_variants, gene_id_to_est_eqtl_effects[gene_id])
			# Borzoi
			borzoi_effects, borzoi_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_est_borzoi_effects[gene_id])
			# Anno
			variant_anno, borzoi_anno_variant_alleles = load_in_snp_gene_data(ordered_cis_variants, gene_id_to_variant_gene_anno[gene_id])

			# Load in LD
			cis_genotype_indices = []
			for var_index, cis_variant in enumerate(ordered_cis_variants):
				cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
				snp_info = rsid_to_snp_info[cis_variant]
				geno_alleles = snp_info[:2]
				
				# Also flip signs of eqtls to match LD
				if np.isnan(eqtl_effects[var_index]) == False:
					if eqtl_variant_alleles[var_index,:][0] == geno_alleles[1]:
						eqtl_effects[var_index] = -1.0*eqtl_effects[var_index]
				if np.isnan(borzoi_effects[var_index]) == False:
					if borzoi_variant_alleles[var_index,:][0] == geno_alleles[1]:
						borzoi_effects[var_index] = -1.0*borzoi_effects[var_index]
				if borzoi_variant_alleles[var_index,:][0] != borzoi_anno_variant_alleles[var_index,:][0]:
					print('annotation alllele assumption erororo')
					pdb.set_trace()
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = G[cis_genotype_indices,:].compute()
			LD = np.corrcoef(geno_mat)

			# Subset LD by missingness
			# A. on eQTL end
			observed_eqtl_indices = np.isnan(eqtl_effects) == False
			eqtl_effects = eqtl_effects[observed_eqtl_indices]
			eqtl_effect_ses = eqtl_effect_ses[observed_eqtl_indices]
			LD = LD[observed_eqtl_indices, :]
			# B. on borzoi end
			observed_borzoi_indices = np.isnan(borzoi_effects) == False
			borzoi_effects = borzoi_effects[observed_borzoi_indices]
			variant_anno = variant_anno[observed_borzoi_indices, :]
			LD = LD[:, observed_borzoi_indices]

			# Add to dictionary
			if gene_id in gene_to_ld_means:
				print('assumption erororo')
				pdb.set_trace()

			gene_to_ld_means[gene_id] = {}
			gene_to_ld_means[gene_id]['eQTL_effect_sizes'] = eqtl_effects
			gene_to_ld_means[gene_id]['eQTL_effect_ses'] = eqtl_effect_ses
			gene_to_ld_means[gene_id]['borzoi_effects'] = borzoi_effects
			gene_to_ld_means[gene_id]['variant_anno'] = variant_anno
			gene_to_ld_means[gene_id]['ld_means'] = LD @ (variant_anno * borzoi_effects[:, None])
			gene_to_ld_means[gene_id]['ld_scores'] = (LD**2) @ variant_anno

	return gene_to_ld_means


def compute_calibration_coefs(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	annos = []
	for gene_name in ordered_gene_names:
		yy.append(gene_id_to_ld_means[gene_name]['eQTL_effect_sizes'])
		xx.append(gene_id_to_ld_means[gene_name]['ld_means'])
		annos.append(gene_id_to_ld_means[gene_name]['variant_anno'])
	yy = np.hstack(yy)
	xx = np.vstack(xx)
	annos = np.vstack(annos)

	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]
	annos = annos[valid_rows, :]


	calibration_coefs, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	
	pred_coefs = np.dot(annos, calibration_coefs)
	avg_pred_coefs = []
	for anno_iter in range(annos.shape[1]):
		averager = np.sum(annos[:, anno_iter]*pred_coefs)/np.sum(annos[:, anno_iter])
		avg_pred_coefs.append(averager)

	return np.asarray(avg_pred_coefs)


def compute_borzoi_variances(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	annos = []
	for gene_name in ordered_gene_names:
		yy.append(np.square(gene_id_to_ld_means[gene_name]['borzoi_effects']))
		xx.append(gene_id_to_ld_means[gene_name]['variant_anno'])
	yy = np.hstack(yy)
	xx = np.vstack(xx)

	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]

	borzoi_var_coefs, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	
	pred_vars = np.dot(xx, borzoi_var_coefs)
	avg_pred_vars = []
	for anno_iter in range(xx.shape[1]):

		averager = np.sum(xx[:, anno_iter]*pred_vars)/np.sum(xx[:, anno_iter])
		avg_pred_vars.append(averager)

	return np.asarray(avg_pred_vars)

def compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names):
	yy = []
	xx = []
	annos = []
	for gene_name in ordered_gene_names:
		yy.append(np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_sizes']) - np.square(gene_id_to_ld_means[gene_name]['eQTL_effect_ses']))
		xx.append(gene_id_to_ld_means[gene_name]['ld_scores'])
		annos.append(gene_id_to_ld_means[gene_name]['variant_anno'])

	yy = np.hstack(yy)
	xx = np.vstack(xx)
	annos = np.vstack(annos)

	valid_rows = np.isfinite(yy) & np.all(np.isfinite(xx), axis=1)
	yy = yy[valid_rows]
	xx = xx[valid_rows, :]
	annos = annos[valid_rows, :]


	taus, _, _, _ = np.linalg.lstsq(xx, yy, rcond=None)
	
	pred_per_snp_h2s = np.dot(annos, taus)
	avg_pred_h2s = []
	for anno_iter in range(annos.shape[1]):
		averager = np.sum(annos[:, anno_iter]*pred_per_snp_h2s)/np.sum(annos[:, anno_iter])
		avg_pred_h2s.append(averager)

	return np.asarray(avg_pred_h2s)


def run_ld_corr(ordered_gene_names, gene_id_to_ld_means):
	avg_calibration_slopes = compute_calibration_coefs(gene_id_to_ld_means, ordered_gene_names)
	avg_borzoi_variances = compute_borzoi_variances(gene_id_to_ld_means, ordered_gene_names)
	avg_per_snp_h2 = compute_eqtl_variances(gene_id_to_ld_means, ordered_gene_names)

	return avg_calibration_slopes, avg_per_snp_h2, avg_borzoi_variances


def create_mapping_from_vg_pair_to_confidently_fine_mapped_eqtl_effect_size(est_fine_mapped_eqtl_effect_size_file, pip_thresh=0.9):
	f = gzip.open(est_fine_mapped_eqtl_effect_size_file, 'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pip = float(data[2])
		if pip < pip_thresh:
			continue
		gene = data[0]
		variant = data[1]
		posterior_mean = float(data[3])
		vg_pair = variant + ':' + gene
		if vg_pair in mapping:
			print('errorr')
			pdb.set_trace()
		mapping[vg_pair] = posterior_mean
	f.close()
	return mapping


def compute_correlation(x, y):
	valid_indices = np.isfinite(x) & np.isfinite(y)
	x = x[valid_indices]
	y = y[valid_indices]
	if len(x) < 2:
		return np.nan
	if np.std(x) == 0.0 or np.std(y) == 0.0:
		return np.nan
	return np.corrcoef(x, y)[0,1]


def compute_regression_slope(x, y):
	valid_indices = np.isfinite(x) & np.isfinite(y)
	x = x[valid_indices]
	y = y[valid_indices]
	if len(x) < 2:
		return np.nan
	x_mean = np.mean(x)
	y_mean = np.mean(y)
	x_var = np.var(x)
	if x_var == 0.0:
		return np.nan
	return np.mean((x - x_mean)*(y - y_mean))/x_var


def compute_fine_mapped_stats(borzoi_effects, eqtl_effects):
	return {
		'correlation': compute_correlation(borzoi_effects, eqtl_effects),
		'calibration_slope': compute_regression_slope(borzoi_effects, eqtl_effects)
	}


def summarize_bootstrap_distribution(observed_value, bootstrap_distribution):
	if np.all(np.isnan(bootstrap_distribution)):
		return np.nan, np.nan, np.nan, np.nan, np.nan
	bootstrap_mean = np.nanmean(bootstrap_distribution)
	bootstrap_se = np.nanstd(bootstrap_distribution, ddof=1)
	if bootstrap_se == 0.0 or np.isnan(bootstrap_se):
		gaussian_z_score = np.nan
	else:
		gaussian_z_score = observed_value/bootstrap_se
	empirical_ci = np.nanquantile(bootstrap_distribution, [.025, .975])
	return bootstrap_mean, bootstrap_se, gaussian_z_score, empirical_ci[0], empirical_ci[1]




##########################
# Command line args
##########################
est_borzoi_effect_size_file = sys.argv[1]
est_fine_mapped_eqtl_effect_size_file = sys.argv[2]
sim_variant_gene_annotation_file = sys.argv[3]
onek_genomes_plink_filestem = sys.argv[4]
fm_corr_output_stem = sys.argv[5]

##############################
# Load in data
##############################

# Create mapping from variant_gene name to confidently fine-mapped eqtl effect sizes
vg_to_fine_mapped_eqtl_effect_size = create_mapping_from_vg_pair_to_confidently_fine_mapped_eqtl_effect_size(est_fine_mapped_eqtl_effect_size_file)

# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file)

# Create mapping from gene id to vector of variant-gene annotations
gene_id_to_variant_gene_anno, anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(sim_variant_gene_annotation_file)


# Initialize vectors to keep track of effects
anno_name_to_borzoi_effect_vec = {}
anno_name_to_eqtl_effect_vec = {}
for anno_name in anno_names:
	anno_name_to_borzoi_effect_vec[anno_name] = []
	anno_name_to_eqtl_effect_vec[anno_name] = []


# Organize
for gene_id in [*gene_id_to_est_borzoi_effects]:
	for var_id in [*gene_id_to_est_borzoi_effects[gene_id]]:
		variant_gene_pair = var_id + ':' + gene_id
		if variant_gene_pair not in vg_to_fine_mapped_eqtl_effect_size:
			continue
		fm_effect_size = vg_to_fine_mapped_eqtl_effect_size[variant_gene_pair]
		borzoi_effect_size = gene_id_to_est_borzoi_effects[gene_id][var_id][6]
		active_annotations = np.where(gene_id_to_variant_gene_anno[gene_id][var_id][6] == 1.0)[0]
		if len(active_annotations) != 1:
			print('annotation assumption error')
			pdb.set_trace()
		anno_name = anno_names[active_annotations[0]]
		anno_name_to_borzoi_effect_vec[anno_name].append(borzoi_effect_size)
		anno_name_to_eqtl_effect_vec[anno_name].append(fm_effect_size)




n_bs = 100
output_names = ['correlation', 'calibration_slope']
observed_stats = {}
bootstrap_stats = {}
for anno_name in anno_names:
	borzoi_effects = np.asarray(anno_name_to_borzoi_effect_vec[anno_name])
	eqtl_effects = np.asarray(anno_name_to_eqtl_effect_vec[anno_name])
	observed_stats[anno_name] = compute_fine_mapped_stats(borzoi_effects, eqtl_effects)
	bootstrap_stats[anno_name] = {}
	for output_name in output_names:
		bootstrap_stats[anno_name][output_name] = np.full(n_bs, np.nan)
	if len(borzoi_effects) < 2:
		continue
	for bs_iter in range(n_bs):
		bs_indices = np.random.choice(np.arange(len(borzoi_effects)), size=len(borzoi_effects), replace=True)
		bs_stats = compute_fine_mapped_stats(borzoi_effects[bs_indices], eqtl_effects[bs_indices])
		for output_name in output_names:
			bootstrap_stats[anno_name][output_name][bs_iter] = bs_stats[output_name]


summary_stats_output_file = fm_corr_output_stem + '_bootstrap_stats.txt'
with open(summary_stats_output_file, 'w') as t:
	t.write('annotation_name\toutput_name\tmean\tbootstrapped_mean\tbootstrap_se\tgaussian_z_score\tempirical_ci_lower\tempirical_ci_upper\n')
	for anno_name in anno_names:
		for output_name in output_names:
			observed_value = observed_stats[anno_name][output_name]
			bootstrap_mean, bootstrap_se, gaussian_z_score, ci_lower, ci_upper = summarize_bootstrap_distribution(observed_value, bootstrap_stats[anno_name][output_name])
			t.write(str(anno_name) + '\t' + output_name + '\t' + str(observed_value) + '\t' + str(bootstrap_mean) + '\t' + str(bootstrap_se) + '\t' + str(gaussian_z_score) + '\t' + str(ci_lower) + '\t' + str(ci_upper) + '\n')

print(summary_stats_output_file)








