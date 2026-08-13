import argparse
import numpy as np
import os
import sys
import pdb
import gzip
from pandas_plink import read_plink


def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')





def create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev, standardize=True):
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

		if var_id not in variant_id_to_genotype_sdev:
			continue

		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()

		geno_sdev = variant_id_to_genotype_sdev[var_id]
		if standardize:
			effect = effect*geno_sdev

		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, effect)
	f.close()
	return mapping


def create_mapping_from_gene_id_to_variant_gene_annotations(sim_variant_gene_annotation_file):
	# New format: columns 6+ are one integer column per annotation, holding the category index
	# the variant-gene pair falls in (-1 if the pair falls in no category of that annotation).
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
		anno = np.asarray(data[6:]).astype(int)
		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('variatn repeat assumption erororo')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, anno)
	f.close()
	return mapping, anno_names


def extract_annotation_categories(annotation_category_file, anno_names):
	# Companion file to the annotation file. Columns (tab-separated, with header):
	# anno_name  source  category_index  category_name
	# Returns, per annotation, the ordered list of its category names.
	anno_name_to_category_names = {}
	anno_name_to_source = {}
	f = open(annotation_category_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_name = data[0]
		source = data[1]
		category_index = int(data[2])
		category_name = data[3]
		if anno_name not in anno_name_to_category_names:
			anno_name_to_category_names[anno_name] = []
			anno_name_to_source[anno_name] = source
		# Categories are written in index order, so appending keeps them aligned to the index
		if category_index != len(anno_name_to_category_names[anno_name]):
			print('assumption eroror: categories not in index order for ' + anno_name)
			pdb.set_trace()
		anno_name_to_category_names[anno_name].append(category_name)
	f.close()

	# Every annotation column in the annotation file needs an entry in the category file
	for anno_name in anno_names:
		if anno_name not in anno_name_to_category_names:
			print('assumption eroror: ' + anno_name + ' missing from category file')
			pdb.set_trace()

	return anno_name_to_category_names, anno_name_to_source


def create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file):	
	variant_id_to_geno_sdev = {}
	f = gzip.open(est_eqtl_effect_size_file,'rt')
	mapping = {}
	head_count = 0
	bad_genes = {}
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
		maf = float(data[9])
		if len(data) < 12 or data[11] == 'nan':
			genotype_sdev = np.sqrt(2.0*maf*(1.0-maf))
		else:
			genotype_sdev = float(data[11])
		std_effect = effect*genotype_sdev
		std_se = se*genotype_sdev
		variant_id_to_geno_sdev[var_id] = genotype_sdev
		if gene_id not in mapping:
			mapping[gene_id] = {}
		if var_id in mapping[gene_id]:
			print('repeat snp error')
			pdb.set_trace()
		mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, pos, a0, a1, std_effect, std_se)
	f.close()
	return mapping, variant_id_to_geno_sdev

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


def load_in_snp_gene_anno_data(ordered_cis_variants, var_to_variant_gene_anno, n_anno):
	# Returns a (n_variant X n_annotation) matrix of category indices. A variant absent from the
	# annotation file gets all -1 (no category in any annotation); such variants are eqtl-only and
	# get subset out on the borzoi end anyways.
	annos = []
	alleles = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_variant_gene_anno:
			annos.append(np.full(n_anno, -1, dtype=int))
			alleles.append(('nan', 'nan'))
		else:
			var_info = var_to_variant_gene_anno[variant_id]
			annos.append(var_info[6])
			alleles.append((var_info[4], var_info[5]))
	return np.vstack(annos), np.asarray(alleles)


def create_one_hot_from_category_indices(category_indices, n_categories):
	# Turn a vector of category indices (-1 = no category) into a one-hot matrix.
	one_hot = np.zeros((len(category_indices), n_categories))
	for var_iter, category_index in enumerate(category_indices):
		if category_index >= 0:
			one_hot[var_iter, category_index] = 1.0
	return one_hot

def extract_gene_chrom_num(var_id_to_est_borzoi_effects):
	var_id = [*var_id_to_est_borzoi_effects][0]
	chrom_num = var_id_to_est_borzoi_effects[var_id][2]
	return chrom_num

def extract_and_write_ld_moments(gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, gene_id_to_variant_gene_anno, genotype_plink_filestem, anno_names, anno_name_to_n_categories, anno_name_to_category_names, genotype_sample_indices, output_file, target_anno_names=None):
	# Stream one line per variant-gene pair with an observed eQTL effect, holding the standardized
	# eQTL effect size alongside, for every annotation category, the LD-mean of the standardized
	# borzoi effects (i.e. the calibration regression's design matrix in sldmc.py, before it gets
	# aggregated to per-gene sufficient statistics).
	# target_anno_names optionally restricts which annotations get LD-mean columns (default: all).
	# The annotation matrix's columns follow anno_names (the annotation file's column order), so
	# each target is looked up by its position in anno_names.
	if target_anno_names is None:
		target_anno_names = anno_names
	target_anno_indices = []
	for target_anno_name in target_anno_names:
		matching_indices = np.where(np.asarray(anno_names) == target_anno_name)[0]
		if len(matching_indices) != 1:
			print('assumption eroror: ' + target_anno_name + ' not found in annotation file')
			pdb.set_trace()
		target_anno_indices.append(matching_indices[0])

	t = gzip.open(output_file, 'wt')
	header_columns = ['gene_id', 'variant_id', 'std_eqtl_effect_size']
	for anno_name in target_anno_names:
		for category_name in anno_name_to_category_names[anno_name]:
			header_columns.append(anno_name + 'X' + category_name)
	t.write('\t'.join(header_columns) + '\n')

	# Loop through chromsomes
	for chrom_num in range(1,23):
		print(chrom_num)

		##################################
		# Load in per-chrom-genotype data
		##################################
		# string of chromosome name
		chrom_string = 'chr' + str(chrom_num)
		# Load in chromosome plink data
		(bim, fam, G) = read_plink(genotype_plink_filestem + str(chrom_num))
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
			variant_anno, borzoi_anno_variant_alleles = load_in_snp_gene_anno_data(ordered_cis_variants, gene_id_to_variant_gene_anno[gene_id], len(anno_names))

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
			
			# Extract genotype
			cis_genotype_indices = np.asarray(cis_genotype_indices)
			# Extract genotype matrix
			geno_mat = (G[cis_genotype_indices,:].compute())[:, genotype_sample_indices]
			row_means = np.nanmean(geno_mat, axis=1)
			nan_rows, nan_cols = np.where(np.isnan(geno_mat))
			geno_mat[nan_rows, nan_cols] = row_means[nan_rows]

			# A variant with no genotype variance in this sample gets an all-nan row AND column in
			# the correlation matrix. The row is harmless (it drops out per-variant below), but a
			# surviving nan column propagates through every matrix product and silently voids the
			# entire gene, so mark these missing on both ends and let the subsetting remove them.
			degenerate_indices = (np.std(geno_mat, axis=1) > 0.0) == False
			eqtl_effects[degenerate_indices] = np.nan
			borzoi_effects[degenerate_indices] = np.nan

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

			# LD-means for every target annotation (columns hstacked in header order)
			ld_mean_columns = []
			for anno_iter, anno_name in zip(target_anno_indices, target_anno_names):
				n_categories = anno_name_to_n_categories[anno_name]
				# One-hot the variants over this annotation's categories (-1 index -> all-zero row)
				one_hot = create_one_hot_from_category_indices(variant_anno[:, anno_iter], n_categories)
				ld_mean_columns.append(LD @ (one_hot * borzoi_effects[:, None]))
			ld_means = np.hstack(ld_mean_columns)

			# Write one line per eQTL-observed variant with finite LD-means in every column
			ordered_eqtl_variants = ordered_cis_variants[observed_eqtl_indices]
			valid_rows = np.isfinite(eqtl_effects) & np.all(np.isfinite(ld_means), axis=1)
			for row_iter in np.where(valid_rows)[0]:
				t.write(gene_id + '\t' + ordered_eqtl_variants[row_iter] + '\t' + str(eqtl_effects[row_iter]) + '\t' + '\t'.join(ld_means[row_iter, :].astype(str)) + '\n')

	t.close()


def create_binned_ld_moment_file(ld_moment_output_file, binned_ld_moment_output_file, n_bins, ld_moment_column_name='interceptXintercept'):
	# Bin the variant-gene pairs from the per-pair LD-moment file into n_bins equally sized groups
	# ordered by their intercept LD-moment, and write one line per bin holding the bin number, the
	# bin's average LD-moment, and the bin's average standardized eQTL effect size.
	ld_moments = []
	eqtl_effects = []
	f = gzip.open(ld_moment_output_file, 'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = np.asarray(data)
			ld_moment_column_indices = np.where(header == ld_moment_column_name)[0]
			eqtl_column_indices = np.where(header == 'std_eqtl_effect_size')[0]
			if len(ld_moment_column_indices) != 1 or len(eqtl_column_indices) != 1:
				print('assumption eroror: required columns missing from ' + ld_moment_output_file)
				pdb.set_trace()
			ld_moment_column_index = ld_moment_column_indices[0]
			eqtl_column_index = eqtl_column_indices[0]
			continue
		ld_moments.append(float(data[ld_moment_column_index]))
		eqtl_effects.append(float(data[eqtl_column_index]))
	f.close()
	ld_moments = np.asarray(ld_moments)
	eqtl_effects = np.asarray(eqtl_effects)

	# Order pairs by LD-moment and split into n_bins equally sized groups (bin sizes differ by at
	# most one pair when the number of pairs is not divisible by n_bins)
	ordering = np.argsort(ld_moments)
	bin_index_groups = np.array_split(ordering, n_bins)

	t = open(binned_ld_moment_output_file, 'w')
	t.write('bin_number\tavg_ld_moment\tavg_std_eqtl_effect_size\n')
	for bin_number, bin_indices in enumerate(bin_index_groups):
		t.write(str(bin_number) + '\t' + str(np.mean(ld_moments[bin_indices])) + '\t' + str(np.mean(eqtl_effects[bin_indices])) + '\n')
	t.close()
	print(binned_ld_moment_output_file)

#########################
# Command line args
#########################
parser = argparse.ArgumentParser(description='Extract LD moments.')
parser.add_argument('--est-borzoi-effect-size-file', dest='est_borzoi_effect_size_file', required=True, help='Estimated Borzoi effect sizes file.')
parser.add_argument('--est-eqtl-effect-size-file', dest='est_eqtl_effect_size_file', required=True, help='Estimated eQTL effect sizes file.')
parser.add_argument('--variant-gene-annotation-file', dest='variant_gene_annotation_file', required=True, help='Variant-gene annotation file.')
parser.add_argument('--genotype-plink-filestem', dest='genotype_plink_filestem', required=True, help='Genotype plink filestem (per-chromosome number appended).')
parser.add_argument('--genotype-sample-mapping-file', dest='genotype_sample_mapping_file', required=True, help='Genotype sample indices for in-sample LD.')
parser.add_argument('--ld-moment-output-stem', dest='ld_moment_output_stem', required=True, help='Output filestem.')
args = parser.parse_args()

est_borzoi_effect_size_file = args.est_borzoi_effect_size_file
est_eqtl_effect_size_file = args.est_eqtl_effect_size_file
variant_gene_annotation_file = args.variant_gene_annotation_file
genotype_plink_filestem = args.genotype_plink_filestem
genotype_sample_mapping_file = args.genotype_sample_mapping_file
ld_moment_output_stem = args.ld_moment_output_stem


# Companion file describing the (annotation, category) pairs in the annotation file
annotation_category_file = variant_gene_annotation_file.split('.txt.gz')[0] + '_categories.txt'

##############################
# Load in data
##############################
# Create mapping from gene id to vector of est (standardized) eqtl effect sizes
gene_id_to_est_eqtl_effects, variant_id_to_genotype_sdev = create_mapping_from_gene_id_to_est_eqtl_effect_sizes(est_eqtl_effect_size_file)


# Create mapping from gene id to vector of est borzoi effects
gene_id_to_est_borzoi_effects = create_mapping_from_gene_id_to_causal_effects(est_borzoi_effect_size_file, variant_id_to_genotype_sdev,standardize=True)


# Create mapping from gene id to vector of variant-gene annotations
gene_id_to_variant_gene_anno, anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(variant_gene_annotation_file)

# Extract the categories making up each annotation
anno_name_to_category_names, anno_name_to_source = extract_annotation_categories(annotation_category_file, anno_names)
anno_name_to_n_categories = {}
for anno_name in anno_names:
	anno_name_to_n_categories[anno_name] = len(anno_name_to_category_names[anno_name])

# Load in genotype sample indices (for this tissue) to achieve in sample ld
genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

##############################
# Extract per variant-gene pair LD-means and write to output file
##############################
# Only the intercept and the (coarse) borzoi magnitude bin annotations get LD-mean columns
target_anno_names = ['intercept', 'borzoi_magnitude_bins']
ld_moment_output_file = ld_moment_output_stem + '_ld_moments.txt.gz'
extract_and_write_ld_moments(gene_id_to_est_borzoi_effects, gene_id_to_est_eqtl_effects, gene_id_to_variant_gene_anno, genotype_plink_filestem, anno_names, anno_name_to_n_categories, anno_name_to_category_names, genotype_sample_indices, ld_moment_output_file, target_anno_names=target_anno_names)


##############################
# Bin variant-gene pairs by their intercept LD-moment and average within bins
##############################
ld_moment_output_file = ld_moment_output_stem + '_ld_moments.txt.gz'
n_bins=100
binned_ld_moment_output_file = ld_moment_output_stem + '_binned_' + str(n_bins) + '_ld_moments.txt'
create_binned_ld_moment_file(ld_moment_output_file, binned_ld_moment_output_file, n_bins)
