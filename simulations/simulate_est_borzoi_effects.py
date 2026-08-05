import numpy as np
import os
import sys
import pdb
import gzip


def sample_borzoi_error_from_simple_mog(pis, sig_sqs):
	component_iter = np.random.choice(np.arange(len(pis)), p=pis)
	error = np.random.normal(loc=0, scale=np.sqrt(sig_sqs[component_iter]))
	return error


def sample_borzoi_effects_from_simple_mog(pis, sig_sqs, causal_variant_gene_effect_size_file, est_borzoi_effect_size_file, sim_borzoi_error_file):
	f = gzip.open(causal_variant_gene_effect_size_file,'rt')
	t = gzip.open(est_borzoi_effect_size_file, 'wt')
	t2 = gzip.open(sim_borzoi_error_file, 'wt')
	t.write('gene\tvariant\tborzoi_effect_size\n')
	t2.write('gene\tvariant\tborzoi_sampling_error\n')

	head_count = 0
	aa = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		variant_id = data[1]
		effect_size = float(data[2])


		borzoi_error = sample_borzoi_error_from_simple_mog(pis, sig_sqs)

		borzoi_effect_size = effect_size + borzoi_error


		t.write(gene_id + '\t' + variant_id + '\t' + str(borzoi_effect_size) + '\n')
		t2.write(gene_id + '\t' + variant_id + '\t' + str(borzoi_error) + '\n')

	f.close()
	t.close()
	t2.close()
	return









#####################
# Command line args
#####################
causal_variant_gene_effect_size_file = sys.argv[1]
est_borzoi_effect_size_file = sys.argv[2]
borzoi_error_distribution = sys.argv[3]
simulation_iter = int(sys.argv[4])
sampling_variance = float(sys.argv[5])
sim_borzoi_error_file = sys.argv[6]


# Set random seed
np.random.seed(simulation_iter)

# Set up sampling variance
if borzoi_error_distribution == 'simple_mog':
	pis = np.asarray([.4, .6])
	sig_sqs = np.asarray([.001, .999])
	scaler = sampling_variance/np.sum(pis*sig_sqs)
	sig_sqs = sig_sqs*scaler

	sample_borzoi_effects_from_simple_mog(pis, sig_sqs, causal_variant_gene_effect_size_file, est_borzoi_effect_size_file, sim_borzoi_error_file)
