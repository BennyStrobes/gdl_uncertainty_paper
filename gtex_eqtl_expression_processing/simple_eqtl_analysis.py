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



def extract_dictionary_list_of_gtf_genes(pc_genes_gtf):
	valid_chroms = {}
	for chrom_num in range(1,23):
		valid_chroms['chr' + str(chrom_num)] = 1


	f = open(pc_genes_gtf)

	dicti = {}
	for line in f:
		if line.startswith('#'):
			continue
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in valid_chroms:
			continue
		if data[2] != 'gene':
			continue

		# Parse attribute field (e.g. 'gene_id "ENSG..."; gene_type "protein_coding"; ...')
		attributes = {}
		for attr in data[8].split(';'):
			attr = attr.strip()
			if attr == '':
				continue
			key, value = attr.split(' ', 1)
			attributes[key] = value.strip('"')

		# Only keep protein coding genes
		'''
		if attributes['gene_type'] != 'protein_coding':
			continue
		'''
			
		ens_id = attributes['gene_id']
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		if ens_id.split('.')[0] in dicti:
			print('repeat gene error')
			pdb.set_trace()
		# TSS is gene start on + strand and gene end on - strand (1-based gtf coordinates)
		strand = data[6]
		if strand == '+':
			tss = int(data[3])
		elif strand == '-':
			tss = int(data[4])
		else:
			print('strand assumption error')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = (data[0], tss, strand)

	f.close()

	return dicti





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

def count_bim_col5_alleles_from_raw_bed(bed_file, n_samples, variant_index):
	# Decode one variant directly from the plink1 .bed (variant-major) file.
	# 2-bit codes (low bits first): 00 = hom bim col5 allele, 01 = missing, 10 = het, 11 = hom bim col6 allele
	bytes_per_variant = (n_samples + 3) // 4
	with open(bed_file, 'rb') as fh:
		magic = fh.read(3)
		if magic != b'\x6c\x1b\x01':
			print('assumption error: unexpected .bed magic bytes / not variant-major')
			pdb.set_trace()
		fh.seek(3 + variant_index*bytes_per_variant)
		raw = np.frombuffer(fh.read(bytes_per_variant), dtype=np.uint8)
	codes = np.stack([(raw >> (2*k)) & 3 for k in range(4)], axis=1).reshape(-1)[:n_samples]
	lookup = np.asarray([2.0, np.nan, 1.0, 0.0])
	return lookup[codes]


def check_alt_allele_dosage_direction(bed_file, G, bim, n_check=20):
	# Verifies that (2.0 - G) equals the count of the ALT allele (from the variant id) for a spread of variants
	n_variants = G.shape[0]
	n_samples = G.shape[1]
	check_indices = np.unique(np.linspace(0, n_variants - 1, n_check).astype(int))
	bim_snps = np.asarray(bim['snp'])
	bim_a0 = np.asarray(bim['a0'])
	for variant_index in check_indices:
		# bim col5 (pandas_plink a0) should be the ALT allele in the variant id
		if bim_a0[variant_index] != bim_snps[variant_index].split('_')[3]:
			print('allele direction error: bim col5 is not ALT for ' + bim_snps[variant_index])
			pdb.set_trace()
		raw_alt_count = count_bim_col5_alleles_from_raw_bed(bed_file, n_samples, variant_index)
		script_alt_count = 2.0 - G[variant_index, :].compute()
		if np.array_equal(raw_alt_count, script_alt_count, equal_nan=True) == False:
			print('allele direction error: 2.0 - G does not equal ALT allele count for ' + bim_snps[variant_index])
			pdb.set_trace()
	return


def extract_snv_mask(bim_snps, bim_a0, bim_a1):
	# True for SNVs (single-base REF and ALT, both in ACGT); False for indels / other variants
	valid_bases = {'A':1, 'C':1, 'G':1, 'T':1}
	is_snv = np.zeros(len(bim_snps), dtype=bool)
	for ii, snp_id in enumerate(bim_snps):
		snp_info = snp_id.split('_')
		ref = snp_info[2]
		alt = snp_info[3]
		if ref != bim_a1[ii] or alt != bim_a0[ii]:
			print('assumption error: variant id alleles do not match bim alleles for ' + snp_id)
			pdb.set_trace()
		if len(ref) == 1 and len(alt) == 1 and ref in valid_bases and alt in valid_bases:
			is_snv[ii] = True
	return is_snv


def extract_ordered_variants_in_cis_window(bim_snps, bim_pos, bim_is_snv, tss, window_size=100000):
	# All SNVs (on this chromosome's plink file) within +/- window_size bp of the TSS, in genotype-file order
	in_window = (np.abs(bim_pos - tss) <= window_size) & bim_is_snv
	return bim_snps[in_window], np.where(in_window)[0]


########################
# Command line args
#########################
expression_file = sys.argv[1]
genotype_sample_mapping_file = sys.argv[2]
plink_stem = sys.argv[3]
gene_gtf_file = sys.argv[4]
new_eqtl_sumstats_output_file = sys.argv[5]


# Open output file handle
t = gzip.open(new_eqtl_sumstats_output_file,'wt')
t.write('gene\tvariant\tchr\tsnp_pos\tA0\tA1\teqtl_effect_size\teqtl_effect_size_se\tN\taf\tdist_to_tss\tgenotype_sdev\n')

# Extract dictionary list of genes
gene_info_dicti = extract_dictionary_list_of_gtf_genes(gene_gtf_file)

# Extract genotype-sample indices corresponding to this tissue
ordered_genotype_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

for chrom_num in range(1,23):

	print(chrom_num, flush=True)

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
	bim_snps = np.asarray(bim['snp'])
	bim_pos = np.asarray(bim['pos']).astype(int)
	# Only test SNVs (filter out indels)
	bim_is_snv = extract_snv_mask(bim_snps, np.asarray(bim['a0']), np.asarray(bim['a1']))
	print('SNVs: ' + str(np.sum(bim_is_snv)) + ' of ' + str(len(bim_is_snv)) + ' variants', flush=True)
	# Confirm that (2.0 - G), used below, counts the ALT allele (= A1 column written to output)
	check_alt_allele_dosage_direction(plink_stem + str(chrom_num) + '.bed', G, bim)


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

		if short_gene_id not in gene_info_dicti:
			print('assumption error: gene not in gtf: ' + short_gene_id)
			pdb.set_trace()
		gene_info = gene_info_dicti[short_gene_id]
		tss = gene_info[1]
		gene_strand = gene_info[2]
		gene_chrom_num = gene_info[0]
		if line_chrom != gene_chrom_num:
			print('assumption error: gene chrom should match line chrom')
			pdb.set_trace()


		# Extract relevent fields for gene
		tss2 = int(data[2])

		expr_vec = np.asarray(data[4:]).astype(float)

		# Extract all SNVs within 100KB of the TSS
		ordered_cis_variants, cis_genotype_indices = extract_ordered_variants_in_cis_window(bim_snps, bim_pos, bim_is_snv, tss, window_size=100000)
		if len(ordered_cis_variants) == 0:
			print('no cis variants for gene: ' + short_gene_id)
			continue

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
