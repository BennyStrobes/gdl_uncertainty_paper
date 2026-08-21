import numpy as np
import os
import sys
import pdb







def load_in_expression_samples(residual_expression_file):
	f = open(residual_expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1

			ordered_samples = np.asarray(data[4:])
			break
	f.close()

	return ordered_samples



###################
# Command line args
###################
residual_expression_file = sys.argv[1]
plink_fam_file = sys.argv[2]
genotype_mapping_file = sys.argv[3]


# Create mapping from gtex-individual ID to genotype index
gtex_id_to_genotype_index = {}
f = open(plink_fam_file)
line_counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')

	gtex_ind_id = data[1].split('-')[0] + '-' + data[1].split('-')[1]

	if gtex_ind_id in gtex_id_to_genotype_index:
		print('assumptioneronrornr')
		pdb.set_trace()
	gtex_id_to_genotype_index[gtex_ind_id] = line_counter
	line_counter = line_counter + 1
f.close()


# Now load in ordered expression samples
tissue_ordered_expresssion_samples = load_in_expression_samples(residual_expression_file)



t = open(genotype_mapping_file,'w')

for sample_id in tissue_ordered_expresssion_samples:

	genotype_index = gtex_id_to_genotype_index[sample_id]
	t.write(str(genotype_index) + '\n')

t.close()

