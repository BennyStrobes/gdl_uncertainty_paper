import numpy as np
import os
import sys
import pdb
import gzip
import glob
import h5py




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

def extract_dictionary_of_borzoi_genes(borzoi_results_dir, chunk_size=5000000):
	# Same borzoi output files as sldmc_analysis/preprocess_borzoi_data_for_sldmc.py
	file_pattern = os.path.join(borzoi_results_dir, 'model_0_chunk_*_borzoi_gtex_only_results.h5')
	borzoi_input_files = sorted(glob.glob(file_pattern))
	if len(borzoi_input_files) == 0:
		print('assumption error: no borzoi result files found matching ' + file_pattern)
		pdb.set_trace()

	dicti = {}
	for borzoi_input_file in borzoi_input_files:
		with h5py.File(borzoi_input_file, 'r') as f:
			gene_ds = f['gene']
			n_rows = gene_ds.shape[0]
			# Only read the gene column, in slices (one row per variant-gene pair, so this is large)
			for start in range(0, n_rows, chunk_size):
				end = min(start + chunk_size, n_rows)
				unique_genes = np.unique(gene_ds[start:end])
				for gene_id in unique_genes:
					if isinstance(gene_id, bytes):
						gene_id = gene_id.decode('utf-8')
					if gene_id.startswith('ENSG') == False:
						print('assumption error: unexpected borzoi gene id ' + gene_id)
						pdb.set_trace()
					dicti[gene_id.split('.')[0]] = 1
	print('number of genes with a borzoi prediction: ' + str(len(dicti)))
	return dicti



#########################
# Command line args
#########################
input_expression_file = sys.argv[1]
input_covariate_file = sys.argv[2]
residual_expression_file = sys.argv[3]
genes_gtf = sys.argv[4]
gtex_donor_info_file = sys.argv[5]
genotype_fam_file = sys.argv[6]
ancestry_label_str = sys.argv[7]
borzoi_results_dir = sys.argv[8]



borzoi_genes_dictionary = extract_dictionary_of_borzoi_genes(borzoi_results_dir)
#protein_coding_genes_dictionary = extract_dictionary_list_of_protein_coding_genes(genes_gtf)
gtex_id_to_ancestry_labels = load_in_ancestry_labels(gtex_donor_info_file)
genotyped_individuals = load_in_genotyped_individuals(genotype_fam_file)

# Extract dictionary list of genes
gene_info_dicti = extract_dictionary_list_of_gtf_genes(genes_gtf)

###############################
# First load in covariate file
cov_mat, ordered_cov_sample_names, ordered_cov_names = load_in_covariate_data(input_covariate_file)

valid_chroms = {}
for chrom_num in range(1,23):
	valid_chroms['chr' + str(chrom_num)] = 1


# Filter to european samples
european_samples = []
for sample_name in ordered_cov_sample_names:
	if gtex_id_to_ancestry_labels[sample_name] == ancestry_label_str:
		if sample_name not in genotyped_individuals:
			european_samples.append(False)
			print('missing genotyped individual: ' + sample_name)
		else:
			european_samples.append(True)
	else:
		european_samples.append(False)
european_samples = np.asarray(european_samples)
print(np.sum(european_samples))
print(np.sum(european_samples)/len(european_samples))

cov_mat = cov_mat[:, european_samples]
ordered_cov_sample_names = ordered_cov_sample_names[european_samples]

cov_design = np.transpose(cov_mat)
intercept = np.ones((cov_design.shape[0], 1))
cov_design = np.hstack((intercept, cov_design))


###############################
# Second, residualize out expression
f = gzip.open(input_expression_file, 'rt')
t = open(residual_expression_file,'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		ordered_expr_sample_names = np.asarray(data[4:])[european_samples]
		if np.array_equal(ordered_expr_sample_names, ordered_cov_sample_names) == False:
			print('assumption eroror')
			pdb.set_trace()
		t.write('\t'.join(data[:4]) + '\t' + '\t'.join(ordered_expr_sample_names) + '\n')
		continue

	gene_id = data[3]

	if gene_id.startswith('ENSG') == False:
		print('assumption erororo')
		pdb.set_trace()
	if len(gene_id.split('.')) != 2:
		print('assumptioneorronror')
		pdb.set_trace()

	if gene_id.split('.')[0] not in borzoi_genes_dictionary:
		continue

	if data[0] not in valid_chroms:
		print('assumption errorr!')
		pdb.set_trace()

	if gene_id.split('.')[0] not in gene_info_dicti:
		print('Gene ' + gene_id + ' not found in gtf file. Gene skipped.')
		continue

	expr_vec = (np.asarray(data[4:]).astype(float))[european_samples]
	beta, _, _, _ = np.linalg.lstsq(cov_design, expr_vec, rcond=None)
	pred_expr = np.dot(cov_design, beta)
	resid_expr = expr_vec - pred_expr

	resid_std = np.std(resid_expr)
	if resid_std == 0.0:
		print('errorororoororororo')
		standardized_resid_expr = np.zeros(len(resid_expr))
	else:
		standardized_resid_expr = (resid_expr - np.mean(resid_expr))/resid_std
	t.write('\t'.join(data[:4]) + '\t' + '\t'.join(standardized_resid_expr.astype(str)) + '\n')
f.close()
t.close()
