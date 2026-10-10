import gzip
import numpy as np
import sys
import pdb

# Computes the true (simulated) correlation and calibration slope between estimated borzoi effects
# and true causal effects, stratified by (annotation, fm_status) pairs, where fm_status is
# 1 if the variant-gene pair is confidently fine-mapped (SuSiE PIP >= pip_thresh) and 0 otherwise
# (including pairs absent from the SuSiE file). Also reports each fm_status pooled across annotations.


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


def load_susie_pip_file(susie_fine_mapping_file):
	mapping = {}
	with open_possibly_gzipped_file(susie_fine_mapping_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				if len(data) < 3 or data[2] != 'susie_pip':
					raise ValueError('Unexpected header in ' + susie_fine_mapping_file + ': ' + '\t'.join(data))
				head_count = head_count + 1
				continue
			key = (data[0], data[1])
			if key in mapping:
				raise ValueError('Duplicate gene-variant pair detected: ' + data[0] + '\t' + data[1])
			mapping[key] = float(data[2])
	return mapping


def load_annotation_file(sim_variant_gene_annotation_file):
	annotation_names = None
	key_to_annotation_name = {}
	with open_possibly_gzipped_file(sim_variant_gene_annotation_file, 'rt') as f:
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				annotation_names = data[6:]
				head_count = head_count + 1
				continue
			key = (data[0], data[1])
			annotation_vec = np.asarray(data[6:]).astype(float)
			active_annotations = np.where(annotation_vec == 1.0)[0]
			if len(active_annotations) != 1:
				raise ValueError('Expected exactly one active annotation for ' + data[0] + '\t' + data[1])
			key_to_annotation_name[key] = annotation_names[active_annotations[0]]
	return annotation_names, key_to_annotation_name


def compute_pearson_correlation(x, y):
	if len(x) < 2:
		return np.nan
	if np.std(x) == 0.0 or np.std(y) == 0.0:
		return np.nan
	return np.corrcoef(x, y)[0,1]


def compute_regression_slope_and_intercept(x, y):
	if len(x) < 2:
		return np.nan, np.nan
	x_mean = np.mean(x)
	y_mean = np.mean(y)
	x_var = np.var(x)
	if x_var == 0.0:
		return np.nan, np.nan
	slope = np.mean((x - x_mean) * (y - y_mean)) / x_var
	intercept = y_mean - slope * x_mean
	return slope, intercept


def write_summary_row(t, annotation_name, fm_status, fm_status_name, est_borzoi, causal):
	est_borzoi = np.asarray(est_borzoi)
	causal = np.asarray(causal)
	pearson_correlation = compute_pearson_correlation(est_borzoi, causal)
	regression_slope, regression_intercept = compute_regression_slope_and_intercept(est_borzoi, causal)
	mean_est_borzoi = np.mean(est_borzoi) if len(est_borzoi) > 0 else np.nan
	mean_causal = np.mean(causal) if len(causal) > 0 else np.nan
	t.write(
		annotation_name + '\t' +
		str(fm_status) + '\t' +
		fm_status_name + '\t' +
		str(len(est_borzoi)) + '\t' +
		str(pearson_correlation) + '\t' +
		str(regression_slope) + '\t' +
		str(regression_intercept) + '\t' +
		str(mean_est_borzoi) + '\t' +
		str(mean_causal) + '\n'
	)


####################
# Command line args
####################
est_borzoi_effect_size_file = sys.argv[1]
sim_variant_gene_annotation_file = sys.argv[2]
causal_variant_gene_effect_size_file = sys.argv[3]
susie_fine_mapping_file = sys.argv[4]
pip_thresh = float(sys.argv[5])
simulation_parameter_summary_file = sys.argv[6]

fm_status_names = ['not_fine_mapped', 'fine_mapped']
pooled_annotation_name = 'all_annotations'


gene_variant_to_est_borzoi = load_effect_size_file(est_borzoi_effect_size_file, 'borzoi_effect_size')
gene_variant_to_causal = load_effect_size_file(causal_variant_gene_effect_size_file, 'effect_size')
gene_variant_to_pip = load_susie_pip_file(susie_fine_mapping_file)
annotation_names, key_to_annotation_name = load_annotation_file(sim_variant_gene_annotation_file)


# Organize effects by (annotation, fm_status)
stratum_to_est_borzoi = {}
stratum_to_causal = {}
for annotation_name in annotation_names + [pooled_annotation_name]:
	for fm_status in [0, 1]:
		stratum_to_est_borzoi[(annotation_name, fm_status)] = []
		stratum_to_causal[(annotation_name, fm_status)] = []

n_pairs = 0
n_fine_mapped = 0
for key in [*key_to_annotation_name]:
	if key not in gene_variant_to_est_borzoi:
		raise ValueError('Missing Borzoi effect size for ' + key[0] + '\t' + key[1])
	if key not in gene_variant_to_causal:
		raise ValueError('Missing causal effect size for ' + key[0] + '\t' + key[1])
	annotation_name = key_to_annotation_name[key]
	# Pairs absent from the SuSiE file are treated as not confidently fine-mapped
	fm_status = 0
	if key in gene_variant_to_pip and gene_variant_to_pip[key] >= pip_thresh:
		fm_status = 1
	n_pairs = n_pairs + 1
	n_fine_mapped = n_fine_mapped + fm_status
	for stratum_annotation_name in [annotation_name, pooled_annotation_name]:
		stratum_to_est_borzoi[(stratum_annotation_name, fm_status)].append(gene_variant_to_est_borzoi[key])
		stratum_to_causal[(stratum_annotation_name, fm_status)].append(gene_variant_to_causal[key])


with open(simulation_parameter_summary_file, 'w') as t:
	t.write('annotation_name\tfm_status\tfm_status_name\tn_variant_gene_pairs\tpearson_correlation\tregression_slope\tregression_intercept\tmean_est_borzoi\tmean_causal\n')
	for annotation_name in annotation_names + [pooled_annotation_name]:
		for fm_status in [0, 1]:
			stratum = (annotation_name, fm_status)
			write_summary_row(t, annotation_name, fm_status, fm_status_names[fm_status], stratum_to_est_borzoi[stratum], stratum_to_causal[stratum])

print(str(n_fine_mapped) + ' of ' + str(n_pairs) + ' variant-gene pairs confidently fine-mapped (PIP >= ' + str(pip_thresh) + ')')
print(simulation_parameter_summary_file)
