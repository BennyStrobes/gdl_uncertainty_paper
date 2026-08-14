import numpy as np
import os
import sys
import pdb
import gzip
import argparse




def extract_expression_difference_result_files(expression_differences_results_dir, tissue_sample_pairs, anno_method):
    """
    File names produced by the difference analysis: one per (unordered) pair of tissues,
    following the naming convention in driver_key.sh.
    """
    result_files = []
    for ii in range(len(tissue_sample_pairs)):
        for jj in range(ii + 1, len(tissue_sample_pairs)):
            target_tissue1, target_sample1 = tissue_sample_pairs[ii].split(':')
            target_tissue2, target_sample2 = tissue_sample_pairs[jj].split(':')
            result_file = expression_differences_results_dir + 'expr_differences_summary_' + target_tissue1 + '_' + target_sample1 + '_vs_' + target_tissue2 + '_' + target_sample2 + '_' + anno_method + '.txt.gz'
            if os.path.exists(result_file) == False:
                print('assumption eroror: missing expression differences result file ' + result_file)
                pdb.set_trace()
            result_files.append(result_file)
    return result_files

def load_in_expression_differences_results(expression_difference_result_files):
    """
    Load in results across all tissue pairs into a single data structure.
    """
    expression_differnce_results = []
    for result_file in expression_difference_result_files:
        with gzip.open(result_file, 'rt') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                data = line.strip().split('\t')
                gene_id = data[0]
                obs_diff = float(data[8])
                pred_diff = float(data[9])
                snr = np.square(pred_diff) / float(data[10])
                expression_differnce_results.append((gene_id, obs_diff, pred_diff, snr))
    return expression_differnce_results


def compute_observed_differences_stratified_by_snr_thresholds(expression_differnce_results, snr_thresholds, n_bootstraps):
    """
    At each SNR threshold, among instances with snr > threshold, compute the average observed
    difference across predicted positives (pred_diff > 0) and across predicted negatives
    (pred_diff < 0). Standard errors come from a gene-level block bootstrap (resample genes
    with replacement; all instances of a resampled gene move together).
    """
    # Organize instance-level data into arrays
    gene_ids = np.asarray([res[0] for res in expression_differnce_results])
    obs_diffs = np.asarray([res[1] for res in expression_differnce_results])
    pred_diffs = np.asarray([res[2] for res in expression_differnce_results])
    snrs = np.asarray([res[3] for res in expression_differnce_results])

    unique_genes, gene_indices = np.unique(gene_ids, return_inverse=True)
    n_genes = len(unique_genes)
    n_thresholds = len(snr_thresholds)

    positive_indices = pred_diffs > 0.0
    negative_indices = pred_diffs < 0.0

    # Per-gene sums and counts of observed differences among predicted positives/negatives, at each threshold
    gene_pos_sums = np.zeros((n_genes, n_thresholds))
    gene_pos_counts = np.zeros((n_genes, n_thresholds))
    gene_neg_sums = np.zeros((n_genes, n_thresholds))
    gene_neg_counts = np.zeros((n_genes, n_thresholds))
    for threshold_iter, snr_threshold in enumerate(snr_thresholds):
        passing_indices = snrs > snr_threshold
        pos = passing_indices & positive_indices
        neg = passing_indices & negative_indices
        gene_pos_sums[:, threshold_iter] = np.bincount(gene_indices[pos], weights=obs_diffs[pos], minlength=n_genes)
        gene_pos_counts[:, threshold_iter] = np.bincount(gene_indices[pos], minlength=n_genes)
        gene_neg_sums[:, threshold_iter] = np.bincount(gene_indices[neg], weights=obs_diffs[neg], minlength=n_genes)
        gene_neg_counts[:, threshold_iter] = np.bincount(gene_indices[neg], minlength=n_genes)

    # Observed averages at each threshold (nan if no instances pass)
    with np.errstate(invalid='ignore', divide='ignore'):
        avg_obs_diff_pos = np.sum(gene_pos_sums, axis=0)/np.sum(gene_pos_counts, axis=0)
        avg_obs_diff_neg = np.sum(gene_neg_sums, axis=0)/np.sum(gene_neg_counts, axis=0)

    # Gene-level block bootstrap
    np.random.seed(0)
    bs_avg_pos = np.zeros((n_bootstraps, n_thresholds))
    bs_avg_neg = np.zeros((n_bootstraps, n_thresholds))
    for bs_iter in range(n_bootstraps):
        # Multiplicity of each gene in this bootstrap sample
        bs_gene_multiplicities = np.bincount(np.random.choice(n_genes, size=n_genes, replace=True), minlength=n_genes)
        with np.errstate(invalid='ignore', divide='ignore'):
            bs_avg_pos[bs_iter, :] = np.dot(bs_gene_multiplicities, gene_pos_sums)/np.dot(bs_gene_multiplicities, gene_pos_counts)
            bs_avg_neg[bs_iter, :] = np.dot(bs_gene_multiplicities, gene_neg_sums)/np.dot(bs_gene_multiplicities, gene_neg_counts)
    avg_obs_diff_pos_se = np.nanstd(bs_avg_pos, axis=0)
    avg_obs_diff_neg_se = np.nanstd(bs_avg_neg, axis=0)

    n_pos_instances = np.sum(gene_pos_counts, axis=0)
    n_neg_instances = np.sum(gene_neg_counts, axis=0)

    return avg_obs_diff_pos, avg_obs_diff_pos_se, avg_obs_diff_neg, avg_obs_diff_neg_se, n_pos_instances, n_neg_instances


parser = argparse.ArgumentParser(description='Merge per-tissue-pair expression differences results.')
parser.add_argument('--expression-differences-results-dir', dest='expression_differences_results_dir', required=True, help='Directory containing the per-tissue-pair expression differences summary files.')
parser.add_argument('--anno-method', dest='anno_method', required=True, help='Name of the annotation method used to generate the borzoi annotation files.')
args = parser.parse_args()

expression_differences_results_dir = args.expression_differences_results_dir
anno_method = args.anno_method

tissue_sample_pairs=["Heart_Left_Ventricle:GTEX-18465-0926-SM-731AY.1", "Brain_Cortex:GTEX-1H3O1-1726-SM-9WYSR.1", "Liver:GTEX-11EQ9-0526-SM-5A5JZ.1", "Whole_Blood:GTEX-1LB8K-0005-SM-DIPED.1", "Muscle_Skeletal:GTEX-13QJ3-0726-SM-5SI68.1"]
n_bootstraps=500


# File names that the difference analysis was run on (one per unordered tissue pair)
expression_difference_result_files = extract_expression_difference_result_files(expression_differences_results_dir, tissue_sample_pairs, anno_method)


# Load in results across all tissue pairs into a single data structure
expression_differnce_results = load_in_expression_differences_results(expression_difference_result_files)


# Continuous grid of SNR thresholds (0 through the 99.99th percentile of observed SNRs)
all_snrs = np.asarray([res[3] for res in expression_differnce_results])
snr_thresholds = np.linspace(0.0, np.quantile(all_snrs, 0.9999), 200)


# Average observed difference among predicted positives/negatives at each SNR threshold (SEs from gene-level block bootstrap)
avg_obs_diff_pos, avg_obs_diff_pos_se, avg_obs_diff_neg, avg_obs_diff_neg_se, n_pos_instances, n_neg_instances = compute_observed_differences_stratified_by_snr_thresholds(expression_differnce_results, snr_thresholds, n_bootstraps)


# Print to output file
output_file = expression_differences_results_dir + 'observed_differences_stratified_by_snr_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('snr_threshold\tavg_observed_diff_predicted_positives\tavg_observed_diff_predicted_positives_se\tavg_observed_diff_predicted_negatives\tavg_observed_diff_predicted_negatives_se\tn_predicted_positives\tn_predicted_negatives\n')
for threshold_iter, snr_threshold in enumerate(snr_thresholds):
    t.write(str(snr_threshold) + '\t' + str(avg_obs_diff_pos[threshold_iter]) + '\t' + str(avg_obs_diff_pos_se[threshold_iter]))
    t.write('\t' + str(avg_obs_diff_neg[threshold_iter]) + '\t' + str(avg_obs_diff_neg_se[threshold_iter]))
    t.write('\t' + str(int(n_pos_instances[threshold_iter])) + '\t' + str(int(n_neg_instances[threshold_iter])) + '\n')
t.close()
print(output_file)
