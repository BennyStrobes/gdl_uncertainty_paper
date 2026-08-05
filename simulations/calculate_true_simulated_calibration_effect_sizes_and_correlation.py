import gzip
import numpy as np
import sys
import pdb

def open_possibly_gzipped_file(file_name, mode):
	if file_name.endswith('.gz'):
		return gzip.open(file_name, mode)
	return open(file_name, mode)


def load_effect_size_file(effect_size_file, effect_size_column_name):
	mapping = {}
	with open_possibly_gzipped_file(effect_size_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				if len(data) < 7 or data[6] != effect_size_column_name:
					raise ValueError('Unexpected header in ' + effect_size_file + ': ' + '\t'.join(data))
				head_count = head_count + 1
				continue
			key = (data[0], data[1])
			if key in mapping:
				raise ValueError('Duplicate gene-variant pair detected: ' + data[0] + '\t' + data[1])
			mapping[key] = float(data[6])
	return mapping


def load_annotation_file(sim_variant_gene_annotation_file):
	annotation_names = None
	annotation_to_keys = {}
	with open_possibly_gzipped_file(sim_variant_gene_annotation_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				annotation_names = data[6:]
				for annotation_name in annotation_names:
					annotation_to_keys[annotation_name] = []
				head_count = head_count + 1
				continue
			key = (data[0], data[1])
			annotation_vec = np.asarray(data[6:]).astype(float)
			active_annotations = np.where(annotation_vec == 1.0)[0]
			if len(active_annotations) != 1:
				raise ValueError('Expected exactly one active annotation for ' + data[0] + '\t' + data[1])
			annotation_name = annotation_names[active_annotations[0]]
			annotation_to_keys[annotation_name].append(key)
	return annotation_names, annotation_to_keys


def compute_regression_slope_and_intercept(x, y):
	x_mean = np.mean(x)
	y_mean = np.mean(y)
	x_var = np.var(x)
	if x_var == 0.0:
		return np.nan, np.nan
	slope = np.mean((x - x_mean) * (y - y_mean)) / x_var
	intercept = y_mean - slope * x_mean
	return slope, intercept


####################
# Command line args
####################
est_borzoi_effect_size_file = sys.argv[1]
sim_variant_gene_annotation_file = sys.argv[2]
causal_variant_gene_effect_size_file = sys.argv[3]
simulation_parameter_summary_file = sys.argv[4]


gene_variant_to_est_borzoi = load_effect_size_file(est_borzoi_effect_size_file, 'borzoi_effect_size')
gene_variant_to_causal = load_effect_size_file(causal_variant_gene_effect_size_file, 'effect_size')
annotation_names, annotation_to_keys = load_annotation_file(sim_variant_gene_annotation_file)


with open(simulation_parameter_summary_file, 'w') as t:
	t.write('annotation_name\tn_variant_gene_pairs\tpearson_correlation\tregression_slope\tregression_intercept\tmean_est_borzoi\tmean_causal\n')
	for annotation_name in annotation_names:
		keys = annotation_to_keys[annotation_name]
		est_borzoi = []
		causal = []
		for key in keys:
			if key not in gene_variant_to_est_borzoi:
				raise ValueError('Missing Borzoi effect size for ' + key[0] + '\t' + key[1])
			if key not in gene_variant_to_causal:
				raise ValueError('Missing causal effect size for ' + key[0] + '\t' + key[1])
			est_borzoi.append(gene_variant_to_est_borzoi[key])
			causal.append(gene_variant_to_causal[key])
		est_borzoi = np.asarray(est_borzoi)
		causal = np.asarray(causal)

		if len(est_borzoi) < 2:
			pearson_correlation = np.nan
		elif np.std(est_borzoi) == 0.0 or np.std(causal) == 0.0:
			pearson_correlation = np.nan
		else:
			pearson_correlation = np.corrcoef(est_borzoi, causal)[0,1]


		regression_slope, regression_intercept = compute_regression_slope_and_intercept(est_borzoi, causal)

		t.write(
			annotation_name + '\t' +
			str(len(est_borzoi)) + '\t' +
			str(pearson_correlation) + '\t' +
			str(regression_slope) + '\t' +
			str(regression_intercept) + '\t' +
			str(np.mean(est_borzoi)) + '\t' +
			str(np.mean(causal)) + '\n'
		)
