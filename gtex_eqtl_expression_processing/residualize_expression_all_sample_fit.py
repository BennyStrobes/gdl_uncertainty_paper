import numpy as np
import os
import sys
import pdb
import gzip


# Part B residualization.
# Per gene, covariate effects are fit on ALL samples of the tissue (all ancestries, using the full
# GTEx v8 covariate set: genotype PCs, all PEER factors, sex, platform, pcr, plus an intercept).
# Residuals are then subset to (a) AFR and (b) EUR genotyped individuals and standardized within each
# ancestry separately. This avoids fitting ~k covariates in a small AFR sample (which attenuates or,
# when n_AFR <= k, destroys the residual signal).
# Genes whose residual sd within an ancestry is ~0 are printed and skipped (not written) for that ancestry.

ZERO_VARIANCE_TOLERANCE = 1e-6


def load_in_covariate_data(input_covariate_file):
	f = open(input_covariate_file)
	sample_names = []
	cov_mat = []
	cov_names = []

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			sample_names = np.asarray(data[1:])
			continue
		cov_names.append(data[0])
		cov_mat.append(np.asarray(data[1:]).astype(float))
	f.close()

	sample_names = np.asarray(sample_names)
	cov_names = np.asarray(cov_names)
	cov_mat = np.asarray(cov_mat)

	return cov_mat, sample_names, cov_names


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
		dicti[ens_id.split('.')[0]] = 1
	f.close()

	return dicti


def load_in_ancestry_labels(gtex_donor_info_file):
	f = open(gtex_donor_info_file)
	mapping = {}
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = np.copy(data)
			continue

		if len(data) != len(header):
			print('assumption eroror')
			pdb.set_trace()
		sample_id = data[0]
		subject_id = sample_id.split('-')[0] + '-' + sample_id.split('-')[1]
		ancestry = data[-1]
		if ancestry == 'NA':
			continue
		if subject_id in mapping:
			if ancestry != mapping[subject_id]:
				print('assumption erororor')
				pdb.set_trace()
		mapping[subject_id] = ancestry
		arr.append(ancestry)
	f.close()
	arr = np.asarray(arr)
	print(np.unique(arr, return_counts=True))
	return mapping


def load_in_genotyped_individuals(genotype_fam_file):
	f = open(genotype_fam_file)
	genotyped_individuals = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		subject_id = data[1]
		if subject_id in genotyped_individuals:
			print('assumption erororor')
			pdb.set_trace()
		genotyped_individuals[subject_id] = 1
	f.close()
	return genotyped_individuals


def get_ancestry_sample_mask(ordered_sample_names, gtex_id_to_ancestry_labels, genotyped_individuals, ancestry_label):
	mask = []
	for sample_name in ordered_sample_names:
		# Samples with NA ancestry are absent from the dictionary; they are used for fitting but never output
		if gtex_id_to_ancestry_labels.get(sample_name, 'NA') == ancestry_label:
			if sample_name not in genotyped_individuals:
				mask.append(False)
				print('missing genotyped individual (' + ancestry_label + '): ' + sample_name)
			else:
				mask.append(True)
		else:
			mask.append(False)
	return np.asarray(mask)


def standardize_or_flag(resid_vec, gene_id, ancestry_label):
	resid_std = np.std(resid_vec)
	if resid_std < ZERO_VARIANCE_TOLERANCE:
		print('zero-variance ' + ancestry_label + ' residual, skipping gene: ' + gene_id + ' (sd=' + str(resid_std) + ')')
		return None
	return (resid_vec - np.mean(resid_vec))/resid_std


#########################
# Command line args
#########################
input_expression_file = sys.argv[1]
input_covariate_file = sys.argv[2]
afr_residual_expression_file = sys.argv[3]
eur_residual_expression_file = sys.argv[4]
gtex_v10_pc_genes_gtf = sys.argv[5]
gtex_donor_info_file = sys.argv[6]
genotype_fam_file = sys.argv[7]

protein_coding_genes_dictionary = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)
gtex_id_to_ancestry_labels = load_in_ancestry_labels(gtex_donor_info_file)
genotyped_individuals = load_in_genotyped_individuals(genotype_fam_file)


###############################
# Load covariates (all samples, full covariate set)
cov_mat, ordered_cov_sample_names, ordered_cov_names = load_in_covariate_data(input_covariate_file)
n_all_samples = len(ordered_cov_sample_names)

cov_design = np.transpose(cov_mat)
intercept = np.ones((cov_design.shape[0], 1))
cov_design = np.hstack((intercept, cov_design))
n_cov_terms = cov_design.shape[1]
cov_rank = np.linalg.matrix_rank(cov_design)

print('Fitting covariates on all ' + str(n_all_samples) + ' samples with ' + str(n_cov_terms) + ' terms (incl. intercept; rank ' + str(cov_rank) + ')')
print('Covariates used: ' + ','.join(ordered_cov_names))
if n_all_samples <= n_cov_terms:
	print('assumption error: number of samples <= number of covariate terms')
	pdb.set_trace()

# Hat-matrix pieces: beta = pinv(C) y (same as least squares solution), computed once for all genes
cov_design_pinv = np.linalg.pinv(cov_design)

# Ancestry masks (over the covariate / expression sample order)
afr_mask = get_ancestry_sample_mask(ordered_cov_sample_names, gtex_id_to_ancestry_labels, genotyped_individuals, 'AFR')
eur_mask = get_ancestry_sample_mask(ordered_cov_sample_names, gtex_id_to_ancestry_labels, genotyped_individuals, 'EUR')
print('AFR samples: ' + str(np.sum(afr_mask)) + ' (' + str(np.sum(afr_mask)/n_all_samples) + ')')
print('EUR samples: ' + str(np.sum(eur_mask)) + ' (' + str(np.sum(eur_mask)/n_all_samples) + ')')
print('Samples with other/NA ancestry used only for fitting: ' + str(n_all_samples - np.sum(afr_mask) - np.sum(eur_mask)))


###############################
# Residualize expression on all samples, then subset + standardize within ancestry
f = gzip.open(input_expression_file, 'rt')
t_afr = open(afr_residual_expression_file, 'w')
t_eur = open(eur_residual_expression_file, 'w')
n_genes_written = 0
n_afr_skipped = 0
n_eur_skipped = 0
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		ordered_expr_sample_names = np.asarray(data[4:])
		if np.array_equal(ordered_expr_sample_names, ordered_cov_sample_names) == False:
			print('assumption eroror: expression and covariate samples differ')
			pdb.set_trace()
		t_afr.write('\t'.join(data[:4]) + '\t' + '\t'.join(ordered_expr_sample_names[afr_mask]) + '\n')
		t_eur.write('\t'.join(data[:4]) + '\t' + '\t'.join(ordered_expr_sample_names[eur_mask]) + '\n')
		continue

	gene_id = data[3]

	if gene_id.startswith('ENSG') == False:
		print('assumption erororo')
		pdb.set_trace()
	if len(gene_id.split('.')) != 2:
		print('assumptioneorronror')
		pdb.set_trace()

	if gene_id.split('.')[0] not in protein_coding_genes_dictionary:
		continue

	expr_vec = np.asarray(data[4:]).astype(float)
	beta = np.dot(cov_design_pinv, expr_vec)
	resid_expr = expr_vec - np.dot(cov_design, beta)

	afr_standardized = standardize_or_flag(resid_expr[afr_mask], gene_id, 'AFR')
	if afr_standardized is None:
		n_afr_skipped = n_afr_skipped + 1
	else:
		t_afr.write('\t'.join(data[:4]) + '\t' + '\t'.join(afr_standardized.astype(str)) + '\n')

	eur_standardized = standardize_or_flag(resid_expr[eur_mask], gene_id, 'EUR')
	if eur_standardized is None:
		n_eur_skipped = n_eur_skipped + 1
	else:
		t_eur.write('\t'.join(data[:4]) + '\t' + '\t'.join(eur_standardized.astype(str)) + '\n')

	n_genes_written = n_genes_written + 1

f.close()
t_afr.close()
t_eur.close()

print('Protein-coding genes processed: ' + str(n_genes_written))
print('Zero-variance genes skipped (AFR): ' + str(n_afr_skipped))
print('Zero-variance genes skipped (EUR): ' + str(n_eur_skipped))
