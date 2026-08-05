import numpy as np
import os
import sys
import pdb
import gzip





def load_in_variants(variant_file):
	f2 = open(variant_file)
	arr = []
	a0_arr = []
	a1_arr = []
	chr_arr = []
	pos_arr = []
	for line in f2:
		line = line.rstrip()
		data = line.split('\t')
		arr.append(data[1])
		a0_arr.append(data[4])
		a1_arr.append(data[5])
		chr_arr.append(data[0])
		pos_arr.append(data[3])
	f2.close()
	return np.asarray(arr), np.asarray(a0_arr), np.asarray(a1_arr), np.asarray(chr_arr), np.asarray(pos_arr)


def simulate_causal_effect_sizes(tot_variants_in_gene, n_causal_variants, cis_snp_h2s):

	if tot_variants_in_gene < 10:
		pdb.set_trace()
	# Initialize causal effeccts 
	causal_effects = np.zeros(tot_variants_in_gene)

	n_gene_causal_variants = np.random.choice(n_causal_variants)

	causal_variant_indices = np.random.choice(np.arange(tot_variants_in_gene), size=n_gene_causal_variants, replace=False)

	for index in causal_variant_indices:
		causal_effects[index] = np.random.normal(loc=0, scale=np.sqrt(np.random.choice(cis_snp_h2s)))


	return causal_effects




########################
# Command line args
########################
simulation_iter = int(sys.argv[1])
gene_ld_summary_file = sys.argv[2]
causal_variant_gene_effect_size_file = sys.argv[3]


# Potential n causal variants per gene
n_causal_variants = np.arange(2,6)
cis_snp_h2s = np.asarray([.015, .02, .025])

# Set random seed
np.random.seed(simulation_iter)

# Open output file handle 
t = gzip.open(causal_variant_gene_effect_size_file,'wt')
t.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\teffect_size\n')


# Loop through genes
f = open(gene_ld_summary_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue

	gene_id = data[0]
	variant_file = data[4]
	cis_variants, cis_variant_a0s, cis_variant_a1s, cis_chroms, cis_poss = load_in_variants(variant_file)

	tot_variants_in_gene = len(cis_variants)

	causal_variant_to_gene_effect_sizes = simulate_causal_effect_sizes(tot_variants_in_gene, n_causal_variants, cis_snp_h2s)
	
	for var_iter, variant_name in enumerate(cis_variants):
		t.write(gene_id + '\t' + variant_name + '\t' + cis_chroms[var_iter] + '\t' + cis_poss[var_iter] + '\t' + cis_variant_a0s[var_iter] + '\t' + cis_variant_a1s[var_iter] + '\t' + str(causal_variant_to_gene_effect_sizes[var_iter]) + '\n')

f.close()
t.close()

