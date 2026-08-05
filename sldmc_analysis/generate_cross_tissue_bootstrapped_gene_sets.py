import numpy as np
import os
import sys
import gzip
import pdb

def extract_tissue_samples_from_target_names_file(borzoi_gtex_independent_target_names_file):
	# Columns (tab-separated, with header):
	# orig_target_index  borzoi_target_index  target_sample  target_description  target_tissue
	tissue_samples = []
	f = open(borzoi_gtex_independent_target_names_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		target_sample = data[2]
		target_tissue = data[4]
		tissue_samples.append((target_tissue, target_sample))
	f.close()
	return tissue_samples


def extract_genes_from_sumstats_file(sumstats_file):
	# Both borzoi effect files and eqtl sumstats files have a header; gene id is in column 0
	genes = {}
	f = gzip.open(sumstats_file, 'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		genes[gene_id] = 1
	f.close()
	return genes


#######################
# Command line args
#######################
eqtl_sumstats_dir = sys.argv[1]
borzoi_output_dir = sys.argv[2]
borzoi_gtex_independent_target_names_file = sys.argv[3]
bootstrapped_cross_tissue_gene_sets_dir = sys.argv[4]


# Parse the (tissue, sample) pairs and build the input filenames for each
tissue_samples = extract_tissue_samples_from_target_names_file(borzoi_gtex_independent_target_names_file)




tissue_to_eqtl_sumstats_file = {}
tissue_to_borzoi_effect_file = {}
tissue_to_genes = {}
for target_tissue, target_sample in tissue_samples:
	print(target_tissue)
	eqtl_sumstats_file = eqtl_sumstats_dir + 'eqtl_results_' + target_tissue + '_sumstats.txt.gz'
	borzoi_effect_file = borzoi_output_dir + target_tissue + '_' + target_sample + '_borzoi_effects.txt.gz'

	if os.path.isfile(borzoi_effect_file) == False:
		print('MISSING borzoi effect file: ' + borzoi_effect_file)
		continue
	if os.path.isfile(eqtl_sumstats_file) == False:
		print('MISSING eqtl sumstats file: ' + eqtl_sumstats_file)
		continue

	tissue_to_eqtl_sumstats_file[target_tissue] = eqtl_sumstats_file
	tissue_to_borzoi_effect_file[target_tissue] = borzoi_effect_file

	# A gene is 'in a tissue' if it appears in BOTH the eqtl sumstats file and the borzoi effect file
	eqtl_genes = set(extract_genes_from_sumstats_file(eqtl_sumstats_file).keys())
	borzoi_genes = set(extract_genes_from_sumstats_file(borzoi_effect_file).keys())
	tissue_to_genes[target_tissue] = eqtl_genes & borzoi_genes

ordered_tissues = sorted(tissue_to_genes.keys())
print(str(len(ordered_tissues)) + ' tissues with input files found')


# Cross-tissue gene set: genes present in at least one tissue
cross_tissue_genes = set()
for target_tissue in ordered_tissues:
	cross_tissue_genes = cross_tissue_genes | tissue_to_genes[target_tissue]

cross_tissue_genes = np.sort(np.asarray(sorted(cross_tissue_genes)))
print(str(len(cross_tissue_genes)) + ' genes present in at least one tissue')


# Write out the (observed) cross-tissue gene set
cross_tissue_gene_set_file = bootstrapped_cross_tissue_gene_sets_dir + 'cross_tissue_gene_set.txt'
t = open(cross_tissue_gene_set_file, 'w')
t.write('gene\n')
for gene_id in cross_tissue_genes:
	t.write(gene_id + '\n')
t.close()


# Generate bootstrapped cross-tissue gene sets (sample genes with replacement)
num_bootstraps = 500
np.random.seed(1)
n_genes = len(cross_tissue_genes)
for bootstrap_iter in range(num_bootstraps):
	sampled_indices = np.random.choice(n_genes, size=n_genes, replace=True)
	bootstrap_gene_set_file = bootstrapped_cross_tissue_gene_sets_dir + 'cross_tissue_gene_set_bootstrap_' + str(bootstrap_iter) + '.txt'
	t = open(bootstrap_gene_set_file, 'w')
	t.write('gene\n')
	for gene_index in sampled_indices:
		t.write(cross_tissue_genes[gene_index] + '\n')
	t.close()

print('done')
