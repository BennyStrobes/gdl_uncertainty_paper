import numpy as np
import os
import sys
import pdb
import argparse


def load_target_names(borzoi_gtex_unique_target_names_file):
	# The borzoi target-names file has a header and five tab-separated columns:
	# orig_target_index  borzoi_target_index  target_identifier  target_description  GTEx_tissue
	# Columns are taken positionally, matching how driver_key.sh reads this same file.
	tissues = []
	samples = []
	head_count = 0
	with open(borzoi_gtex_unique_target_names_file) as f:
		for line in f:
			line = line.rstrip('\n')
			if line == '':
				continue
			data = line.split('\t')
			if len(data) < 5:
				raise ValueError('Expected at least 5 tab-separated columns in ' + borzoi_gtex_unique_target_names_file + ', found ' + str(len(data)))
			if head_count == 0:
				head_count = head_count + 1
				continue
			samples.append(data[2])
			tissues.append(data[4])

	if len(tissues) < 2:
		raise ValueError('Need at least 2 tissues to permute, found ' + str(len(tissues)) + ' in ' + borzoi_gtex_unique_target_names_file)
	# Output stems downstream are keyed on tissue name, so a repeated tissue would collide
	if len(set(tissues)) != len(tissues):
		raise ValueError('Repeated GTEx tissue in ' + borzoi_gtex_unique_target_names_file)

	return tissues, samples


def draw_permuted_partner_indices(n_tissues, seed):
	# For each tissue, draw a partner tissue uniformly at random from the other tissues. Partners are
	# drawn independently, so one tissue can donate its borzoi predictions to several others.
	rng = np.random.default_rng(seed)
	partner_indices = []
	for tissue_iter in range(n_tissues):
		candidate_indices = [other_iter for other_iter in range(n_tissues) if other_iter != tissue_iter]
		partner_indices.append(int(rng.choice(candidate_indices)))
	return partner_indices


#######################
# Command line args
#######################
parser = argparse.ArgumentParser()
parser.add_argument('--borzoi-gtex-unique-target-names-file', dest='borzoi_gtex_unique_target_names_file', required=True, type=str, help='Borzoi gtex target-names file (one row per tissue)')
parser.add_argument('--tissue-permuted-pairs-output-file', dest='tissue_permuted_pairs_output_file', required=True, type=str, help='Output file: one row per tissue, holding that tissue and a randomly drawn partner tissue/sample')
parser.add_argument('--seed', dest='seed', default=0, type=int, help='Random seed, so the drawn pairing is reproducible')
args = parser.parse_args()

borzoi_gtex_unique_target_names_file = args.borzoi_gtex_unique_target_names_file
tissue_permuted_pairs_output_file = args.tissue_permuted_pairs_output_file
seed = args.seed


##############################
# Load tissues and draw a permuted partner for each
##############################
tissues, samples = load_target_names(borzoi_gtex_unique_target_names_file)
print('Loaded ' + str(len(tissues)) + ' tissues from ' + borzoi_gtex_unique_target_names_file)

partner_indices = draw_permuted_partner_indices(len(tissues), seed)


##############################
# Print permuted pairing
##############################
# GTEx_tissue / target_identifier name the tissue supplying the eQTL sumstats and genotype sample
# mapping; permuted_GTEx_tissue / permuted_target_identifier name the tissue supplying the borzoi
# effects and borzoi annotations.
with open(tissue_permuted_pairs_output_file, 'w') as t:
	t.write('GTEx_tissue\ttarget_identifier\tpermuted_GTEx_tissue\tpermuted_target_identifier\n')
	for tissue_iter in range(len(tissues)):
		partner_iter = partner_indices[tissue_iter]
		if partner_iter == tissue_iter:
			raise ValueError('Drew a self-pairing for ' + tissues[tissue_iter])
		t.write(tissues[tissue_iter] + '\t' + samples[tissue_iter] + '\t' + tissues[partner_iter] + '\t' + samples[partner_iter] + '\n')

print(tissue_permuted_pairs_output_file)
