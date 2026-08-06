import numpy as np
import os
import sys
import pdb
import gzip
















#####################
# Command line args
#####################
causal_variant_gene_effect_size_file = sys.argv[1]
est_borzoi_effect_size_file = sys.argv[2]
simulation_iter = int(sys.argv[3])
sim_variant_gene_annotation_file = sys.argv[4]
n_anno = int(sys.argv[5])
sldmc_variant_gene_annotation_file = sys.argv[6]
sldmc_annotation_category_file = sys.argv[7]

np.random.seed(simulation_iter)


f = gzip.open(causal_variant_gene_effect_size_file,'rt')
t = gzip.open(est_borzoi_effect_size_file, 'wt')
t2 = gzip.open(sim_variant_gene_annotation_file, 'wt')
t.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\tborzoi_effect_size\n')
t2.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1')
for ii in range(n_anno):
	t2.write('\t' + 'anno' + str(ii))
t2.write('\n')

gaussian_prop_grid = np.linspace(0, 1, n_anno)

# Same annotations, written in the format the standalone SLDMC package expects: one column per
# annotation holding the 0-based category index, rather than a one-hot column per category.
t3 = gzip.open(sldmc_variant_gene_annotation_file, 'wt')
t3.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\tintercept\tsim_signal_prop\n')

# SLDMC locates this companion file itself, by swapping '.txt.gz' on the annotation file for
# '_categories.txt', so the two names passed in have to follow that convention.
t4 = open(sldmc_annotation_category_file,'w')
t4.write('anno_name\tsource\tcategory_index\tcategory_name\n')
t4.write('intercept\tsimulation\t0\tintercept\n')
for ii in range(n_anno):
	t4.write('sim_signal_prop\tsimulation\t' + str(ii) + '\tprop_' + str(round(gaussian_prop_grid[ii],3)) + '\n')
t4.close()


head_count = 0

gene_to_probs = {}

aa = []
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene_id = data[0]
	variant_id = data[1]
	chromer = data[2]
	snp_pos = data[3]
	a0 = data[4]
	a1 = data[5]
	effect_size = float(data[6])

	anno_vec = np.zeros(n_anno)

	if gene_id not in gene_to_probs:
		probs = np.random.dirichlet(np.ones(n_anno)*.25, size=None)
		gene_to_probs[gene_id] = probs

	anno_index = np.random.choice(n_anno, p=gene_to_probs[gene_id])

	anno_vec[anno_index] = 1


	constant_small_gaussian_error = np.random.normal(loc=0,scale=np.sqrt(.00005))

	gaussian_error = np.random.normal(loc=0,scale=np.sqrt(.005))

	standardized_borzoi_effect = (3.0*gaussian_prop_grid[anno_index]*effect_size + (1.0-gaussian_prop_grid[anno_index])*gaussian_error) + constant_small_gaussian_error

	t.write(gene_id + '\t' + variant_id + '\t' + chromer + '\t' + snp_pos + '\t' + a0 + '\t' + a1 + '\t' + str(standardized_borzoi_effect) + '\n')

	t2.write(gene_id + '\t' + variant_id + '\t' + chromer + '\t' + snp_pos + '\t' + a0 + '\t' + a1 + '\t' + '\t'.join(anno_vec.astype(str)) + '\n')

	t3.write(gene_id + '\t' + variant_id + '\t' + chromer + '\t' + snp_pos + '\t' + a0 + '\t' + a1 + '\t' + '0' + '\t' + str(anno_index) + '\n')


f.close()
t.close()
t2.close()
t3.close()


