import numpy as np
import os
import sys
import pdb
import pickle
import glob
from pandas_plink import read_plink
import pyarrow.parquet as pq
import pandas as pd
import pyarrow.compute as pc
import h5py
import gzip



def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf):
	valid_chroms = {}
	for chrom_num in range(1,23):
		valid_chroms['chr' + str(chrom_num)] = 1


	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in valid_chroms:
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		strand = data[6]
		if strand not in ['+', '-']:
			print('assumption oerororo')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = strand

	f.close()

	return dicti




def load_in_eqtl_data(orig_eqtl_sumstats_h5_stem, chrom_num, pc_genes):
	filer = orig_eqtl_sumstats_h5_stem + str(chrom_num) + '.parquet'
	pf = pq.ParquetFile(filer)
	mapping = {}
	for rg in range(pf.num_row_groups):
		table = pf.read_row_group(
			rg,
			columns=['gene_id', 'variant_id', 'tss_distance', 'af', 'slope', 'slope_se', 'ma_count']
		)

		if table.num_columns == 0:
			continue

		gene_col = table['gene_id']
		variant_col = table['variant_id']
		dist_col = table['tss_distance']
		af_col = table['af']
		slope_col = table['slope']
		slope_se_col = table['slope_se']


		gene_ids = table['gene_id'].to_pylist()
		variant_ids = table['variant_id'].to_pylist()
		tss_dists = table['tss_distance'].to_pylist()
		afs = table['af'].to_pylist()
		ma_counts = table['ma_count'].to_pylist()
		slopes = table['slope'].to_pylist()
		slope_ses = table['slope_se'].to_pylist()

		for ii, gene_id_full in enumerate(gene_ids):
			gene_id = gene_id_full.split('.')[0]
			if gene_id not in pc_genes:
				continue

			maf = afs[ii]
			if maf > 0.5:
				maf = 1.0 - maf
			line_tss_dist = tss_dists[ii]
			if np.abs(line_tss_dist) > 100000:
				continue
			variant_id = variant_ids[ii]
			var_info = variant_id.split('_')
			if len(var_info[2]) != 1 or len(var_info[3]) != 1:
				continue
			if slopes[ii] is None or slope_ses[ii] is None:
				continue
			if slope_ses[ii] == 0:
				continue

			if gene_id not in mapping:
				mapping[gene_id] = {}
			zed = slopes[ii]/slope_ses[ii]
			if variant_id in mapping[gene_id]:
				print('repeat snp error')
				pdb.set_trace()
			mapping[gene_id][variant_id] = (gene_id, variant_id, chrom_num, var_info[1], var_info[2], var_info[3], slopes[ii], slope_ses[ii], afs[ii], line_tss_dist)
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

		snp_info = snp_id.split('_')
		if snp_info[2] != a1_arr[ii]:
			print('error')
			pdb.set_trace()
		if snp_info[3] != a0_arr[ii]:
			print('errror')
			pdb.set_trace()

		dicti[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
	return dicti

def extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_est_eqtl_effects):
	unique_vars = np.unique([*var_to_est_eqtl_effects])
	final_vars = []
	for var in unique_vars:
		if var not in rsid_to_genotype_index:
			continue
		geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])

		passing = False
		if var in var_to_est_eqtl_effects:
			passing = True
			eqtl_alleles = var_to_est_eqtl_effects[var][4:6]
			if set(geno_alleles) != set(eqtl_alleles):
				passing = False

		if passing == False:
			print('hit!!!!')
			continue
		final_vars.append(var)
	return np.asarray(final_vars)


def load_in_snp_gene_eqtl_data(ordered_cis_variants, var_to_est_eqtl_effects):
	effects = []
	alleles = []
	effect_ses = []
	mafs = []
	distys = []

	for variant_id in ordered_cis_variants:
		if variant_id not in var_to_est_eqtl_effects:
			print('assumpriontonoernorr')
			pdb.set_trace()

		var_info = var_to_est_eqtl_effects[variant_id]
		effects.append(var_info[6])
		effect_ses.append(var_info[7])
		alleles.append(var_info[4] + '_' + var_info[5])
		mafs.append(var_info[8])
		distys.append(var_info[9])
	return np.asarray(effects), np.asarray(alleles), np.asarray(effect_ses), np.asarray(mafs), np.asarray(distys)


########################
# Command line args
#########################
expression_file = sys.argv[1]
genotype_sample_mapping_file = sys.argv[2]
plink_stem = sys.argv[3]
orig_eqtl_sumstats_h5_stem = sys.argv[4]
new_eqtl_sumstats_output_file = sys.argv[5]
gtex_v10_pc_genes_gtf = sys.argv[6]


# Open output file handle
t = gzip.open(new_eqtl_sumstats_output_file,'wt')
t.write('gene\tvariant\tchr\tsnp_pos\tA0\tA1\teqtl_effect_size\teqtl_effect_size_se\tN\taf\tdist_to_tss\tgenotype_sdev\n')

# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)

# Extract genotype-sample indices corresponding to this tissue
ordered_genotype_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

for chrom_num in range(1,23):

	print(chrom_num, flush=True)

	####################################
	# First, load in existing eQTL data for this tissue (used to control which vg-pairs we test)
	eqtl_data_obj = load_in_eqtl_data(orig_eqtl_sumstats_h5_stem, chrom_num, pc_genes)

	####################################
	# Second, Load in per-chrom-genotype data
	# string of chromosome name
	chrom_string = 'chr' + str(chrom_num)
	# Load in chromosome plink data
	(bim, fam, G) = read_plink(plink_stem + str(chrom_num))
	# Create mapping from variant id to index
	rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
	# Create mapping from rsid to a0, a1
	rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))


	####################################
	# Third, loop through genes

	f = open(expression_file)
	var_counter = 0
	skip_counter = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Filter out genes not on this chromosome
		line_chrom = data[0]
		if line_chrom != chrom_string:
			continue

		# Filter out genes not in existing eqtl data
		full_gene_id = data[3]
		short_gene_id = full_gene_id.split('.')[0]
		if short_gene_id not in eqtl_data_obj:
			continue

		gene_strand = pc_genes[short_gene_id]

		# Extract relevent fields for gene
		tss = int(data[2])
		expr_vec = np.asarray(data[4:]).astype(float)


		# Extract ordered list of variants used previously for eQTL calling
		ordered_cis_variants = extract_ordered_variants_to_test_on_gene(rsid_to_genotype_index, rsid_to_snp_info, eqtl_data_obj[short_gene_id])


		# Load in existing eqtl data for the gene
		orig_eqtl_effects, orig_eqtl_variant_alleles, orig_eqtl_effect_ses, orig_eqtl_afs, orig_eqtl_distance = load_in_snp_gene_eqtl_data(ordered_cis_variants, eqtl_data_obj[short_gene_id])
		
		cis_genotype_indices = []
		for var_index, cis_variant in enumerate(ordered_cis_variants):
			cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
			snp_info = rsid_to_snp_info[cis_variant]
			geno_alleles = snp_info[:2]
			geno_allele_string = geno_alleles[0] + '_' + geno_alleles[1]
			geno_allele_alt_string = geno_alleles[1] + '_' + geno_alleles[0]

			if geno_allele_alt_string != orig_eqtl_variant_alleles[var_index]:
				print('assumptionerororo')
				pdb.set_trace()


		# Extract genotype
		cis_genotype_indices = np.asarray(cis_genotype_indices)
		# Extract genotype matrix
		geno_mat = (G[cis_genotype_indices,:].compute())[:, ordered_genotype_indices]

		# Get genotype mat
		per_snp_missing_rate = np.sum(np.isnan(geno_mat),axis=1)/(geno_mat.shape[1])

		geno_mat = 2.0 - geno_mat


		ests= []
		ses = []
		afs = []
		# Now loop through variants
		for var_index, cis_variant in enumerate(ordered_cis_variants):

			var_counter = var_counter + 1
			# Throw out snps with any missingness
			if per_snp_missing_rate[var_index] > 0.0:
				skip_counter = skip_counter + 1
				continue

			# Extract distance to TSS
			var_pos = int(cis_variant.split('_')[1])
			dist_to_tss = var_pos - tss

			if gene_strand == '+':
				dist_to_tss = dist_to_tss*1.0
			elif gene_strand == '-':
				dist_to_tss = dist_to_tss*-1.0
			else:
				print('asssumption erroror')
				pdb.set_trace()

			# Extract variant af
			var_af = np.mean(geno_mat[var_index,:])/2.0

			# Compute marginal eQTL effect size and standard error with an intercept
			geno_vec = geno_mat[var_index,:]
			geno_mean = np.mean(geno_vec)
			centered_geno_vec = geno_vec - geno_mean
			sxx = np.sum(np.square(centered_geno_vec))

			N_samp = len(geno_vec)

			snp_info = cis_variant.split('_')
			snp_pos = snp_info[1]
			snp_a0 = snp_info[2]
			snp_a1 = snp_info[3]

			genotype_sdev = np.std(geno_vec)

			if sxx == 0.0 or len(expr_vec) <= 2:
				eqtl_effect = np.nan
				eqtl_effect_se = np.nan
			else:
				eqtl_effect = np.sum(centered_geno_vec*expr_vec)/sxx
				expr_intercept = np.mean(expr_vec) - (eqtl_effect*geno_mean)
				resid_expr = expr_vec - (expr_intercept + (eqtl_effect*geno_vec))
				resid_var = np.sum(np.square(resid_expr))/(len(expr_vec) - 2.0)
				eqtl_effect_se = np.sqrt(resid_var/sxx)


				t.write(short_gene_id + '\t' + cis_variant + '\t' + str(chrom_num) + '\t' + snp_pos + '\t' + snp_a0 + '\t' + snp_a1 + '\t' + str(eqtl_effect) + '\t' + str(eqtl_effect_se) + '\t')
				t.write(str(N_samp) + '\t' + str(var_af) + '\t' + str(dist_to_tss) + '\t' + str(genotype_sdev) + '\n')

	f.close()
	print(skip_counter, var_counter, flush=True)



t.close()
