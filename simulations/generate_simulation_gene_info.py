import numpy as np
import os
import sys
import pdb



def extract_ordered_list_of_genes(gene_annotation_file, chrom_num):
	f = open(gene_annotation_file)
	gene_info_tupler = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		if data[0] != 'chr' + str(chrom_num):
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		strand = data[6]
		if strand == '+':
			tss = data[3]
		else:
			tss = data[4]

		gene_info_tupler.append((ens_id, 'chr' + str(chrom_num), int(tss)))
	f.close()

	return gene_info_tupler

def load_in_variant_info(chrom_bim_file):
	f = open(chrom_bim_file)
	var_info = []
	var_tss = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		var_info.append(np.asarray(data))
		var_tss.append(int(data[3]))
	f.close()

	return np.asarray(var_info), np.asarray(var_tss)




def create_cis_var_bim_file(cis_var_bim_file, cis_var_info):
	t2 = open(cis_var_bim_file,'w')

	for row_iter in range(cis_var_info.shape[0]):
		t2.write('\t'.join(cis_var_info[row_iter,:]) + '\n')
	t2.close()
	return

####################
# Command line args
####################
gene_annotation_file = sys.argv[1]
onek_genomes_plink_filestem = sys.argv[2]
gene_info_dir = sys.argv[3]


# Open gene summary output file
gene_summary_output_file = gene_info_dir + 'gene_summary.txt'
t = open(gene_summary_output_file,'w')
t.write('gene_name\tchrom_num\tgene_tss\tn_cis_variants\tcis_variant_bim_file\n')

# Loop through chromosomes
for chrom_num in range(1,23):
	# Extract ordered list of genes
	gene_tupler = extract_ordered_list_of_genes(gene_annotation_file, chrom_num)

	# Load in variant info
	variant_names, variant_positions = load_in_variant_info(onek_genomes_plink_filestem + str(chrom_num) + '.bim')


	# loop through genes on this chromosome
	for gene_info in gene_tupler:
		gene_name = gene_info[0]
		gene_chrom = gene_info[1]
		gene_tss = gene_info[2]

		cis_var_indices = np.abs(variant_positions - gene_tss) < 100000
		n_cis_variants = np.sum(cis_var_indices)

		# Downstream inference requires a minimum number of cis variants per gene
		if n_cis_variants < 10:
			continue

		cis_var_info = variant_names[cis_var_indices, :]

		cis_var_bim_file = gene_info_dir + gene_name + '_cis_variant_bim_file.txt'
		create_cis_var_bim_file(cis_var_bim_file, cis_var_info)

		t.write(gene_name + '\t' + 'chr' + str(chrom_num) + '\t' + str(gene_tss) + '\t' + str(n_cis_variants) + '\t' + cis_var_bim_file + '\n')

t.close()

