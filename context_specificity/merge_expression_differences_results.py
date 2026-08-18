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
    tissue_pair_labels = []
    for ii in range(len(tissue_sample_pairs)):
        for jj in range(ii + 1, len(tissue_sample_pairs)):
            target_tissue1, target_sample1 = tissue_sample_pairs[ii].split(':')
            target_tissue2, target_sample2 = tissue_sample_pairs[jj].split(':')
            result_file = expression_differences_results_dir + 'expr_differences_summary_' + target_tissue1 + '_' + target_sample1 + '_vs_' + target_tissue2 + '_' + target_sample2 + '_' + anno_method + '.txt.gz'
            if os.path.exists(result_file) == False:
                print('assumption eroror: missing expression differences result file ' + result_file)
                pdb.set_trace()
            result_files.append(result_file)
            tissue_pair_labels.append(target_tissue1 + '_vs_' + target_tissue2)
    return result_files, tissue_pair_labels

def load_in_expression_differences_results(expression_difference_result_files, tissue_pair_labels):
    """
    Load in results across all tissue pairs into a single data structure.
    """
    expression_differnce_results = []
    for result_file, tissue_pair_label in zip(expression_difference_result_files, tissue_pair_labels):
        with gzip.open(result_file, 'rt') as f:
            header = f.readline().strip().split('\t')
            alt_tau_var_index = header.index('predicted_expr_difference_uncertainty_var_alt_tau')
            for line in f:
                data = line.strip().split('\t')
                gene_id = data[0]
                obs_diff = float(data[8])
                pred_diff = float(data[9])
                snr = np.square(pred_diff) / float(data[10])
                # Alternative-tau SNR; a non-positive uncertainty variance (every residual
                # variance truncated to zero) leaves the SNR undefined -> nan (excluded by
                # every threshold mask)
                alt_tau_var = float(data[alt_tau_var_index])
                snr_alt_tau = np.square(pred_diff) / alt_tau_var if alt_tau_var > 0.0 else np.nan
                expression_differnce_results.append((gene_id, obs_diff, pred_diff, snr, snr_alt_tau, tissue_pair_label))
    return expression_differnce_results


def organize_expression_difference_results_into_arrays(expression_differnce_results):
    # Organize instance-level data into arrays
    gene_ids = np.asarray([res[0] for res in expression_differnce_results])
    obs_diffs = np.asarray([res[1] for res in expression_differnce_results])
    pred_diffs = np.asarray([res[2] for res in expression_differnce_results])
    snrs = np.asarray([res[3] for res in expression_differnce_results])
    snr_alt_taus = np.asarray([res[4] for res in expression_differnce_results])
    instance_tissue_pair_labels = np.asarray([res[5] for res in expression_differnce_results])
    return gene_ids, obs_diffs, pred_diffs, snrs, snr_alt_taus, instance_tissue_pair_labels


def compute_observed_differences_over_instance_masks(gene_ids, obs_diffs, pred_diffs, instance_masks, n_bootstraps):
    """
    For each instance mask, among masked instances, compute the average observed difference
    across predicted positives (pred_diff > 0) and across predicted negatives (pred_diff < 0).
    Standard errors come from a gene-level block bootstrap (resample genes with replacement;
    all instances of a resampled gene move together).
    """
    unique_genes, gene_indices = np.unique(gene_ids, return_inverse=True)
    n_genes = len(unique_genes)
    n_masks = len(instance_masks)

    positive_indices = pred_diffs > 0.0
    negative_indices = pred_diffs < 0.0

    # Per-gene sums and counts of observed differences among predicted positives/negatives, per mask
    gene_pos_sums = np.zeros((n_genes, n_masks))
    gene_pos_counts = np.zeros((n_genes, n_masks))
    gene_neg_sums = np.zeros((n_genes, n_masks))
    gene_neg_counts = np.zeros((n_genes, n_masks))
    for mask_iter, instance_mask in enumerate(instance_masks):
        pos = instance_mask & positive_indices
        neg = instance_mask & negative_indices
        gene_pos_sums[:, mask_iter] = np.bincount(gene_indices[pos], weights=obs_diffs[pos], minlength=n_genes)
        gene_pos_counts[:, mask_iter] = np.bincount(gene_indices[pos], minlength=n_genes)
        gene_neg_sums[:, mask_iter] = np.bincount(gene_indices[neg], weights=obs_diffs[neg], minlength=n_genes)
        gene_neg_counts[:, mask_iter] = np.bincount(gene_indices[neg], minlength=n_genes)

    # Observed averages per mask (nan if no instances)
    with np.errstate(invalid='ignore', divide='ignore'):
        avg_obs_diff_pos = np.sum(gene_pos_sums, axis=0)/np.sum(gene_pos_counts, axis=0)
        avg_obs_diff_neg = np.sum(gene_neg_sums, axis=0)/np.sum(gene_neg_counts, axis=0)

    # Gene-level block bootstrap
    np.random.seed(0)
    bs_avg_pos = np.zeros((n_bootstraps, n_masks))
    bs_avg_neg = np.zeros((n_bootstraps, n_masks))
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


def compute_observed_differences_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snrs, snr_thresholds, n_bootstraps):
    """
    At each threshold, among instances whose stratification value (the snrs argument; any
    quantity works) exceeds the threshold, compute the average observed difference across
    predicted positives/negatives.
    """
    instance_masks = [snrs > snr_threshold for snr_threshold in snr_thresholds]
    return compute_observed_differences_over_instance_masks(gene_ids, obs_diffs, pred_diffs, instance_masks, n_bootstraps)


def compute_sign_concordance_over_instance_masks(gene_ids, obs_diffs, pred_diffs, instance_masks, n_bootstraps):
    """
    For each instance mask, among masked instances, compute the fraction of instances where the
    sign of the observed difference matches the sign of the predicted difference.
    Standard errors come from a gene-level block bootstrap.
    """
    unique_genes, gene_indices = np.unique(gene_ids, return_inverse=True)
    n_genes = len(unique_genes)
    n_masks = len(instance_masks)

    concordant_indices = np.sign(obs_diffs) == np.sign(pred_diffs)

    # Per-gene concordant counts and total counts per mask
    gene_concordant_counts = np.zeros((n_genes, n_masks))
    gene_counts = np.zeros((n_genes, n_masks))
    for mask_iter, instance_mask in enumerate(instance_masks):
        gene_concordant_counts[:, mask_iter] = np.bincount(gene_indices[instance_mask & concordant_indices], minlength=n_genes)
        gene_counts[:, mask_iter] = np.bincount(gene_indices[instance_mask], minlength=n_genes)

    with np.errstate(invalid='ignore', divide='ignore'):
        sign_concordances = np.sum(gene_concordant_counts, axis=0)/np.sum(gene_counts, axis=0)

    # Gene-level block bootstrap
    np.random.seed(0)
    bs_concordances = np.zeros((n_bootstraps, n_masks))
    for bs_iter in range(n_bootstraps):
        bs_gene_multiplicities = np.bincount(np.random.choice(n_genes, size=n_genes, replace=True), minlength=n_genes)
        with np.errstate(invalid='ignore', divide='ignore'):
            bs_concordances[bs_iter, :] = np.dot(bs_gene_multiplicities, gene_concordant_counts)/np.dot(bs_gene_multiplicities, gene_counts)
    sign_concordance_ses = np.nanstd(bs_concordances, axis=0)

    n_instances = np.sum(gene_counts, axis=0)

    return sign_concordances, sign_concordance_ses, n_instances


def compute_sign_concordance_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snrs, snr_thresholds, n_bootstraps):
    """
    At each threshold, among instances whose stratification value (the snrs argument; any
    quantity works) exceeds the threshold, compute the sign concordance.
    """
    instance_masks = [snrs > snr_threshold for snr_threshold in snr_thresholds]
    return compute_sign_concordance_over_instance_masks(gene_ids, obs_diffs, pred_diffs, instance_masks, n_bootstraps)


def compute_obs_vs_pred_slope_over_instance_masks(gene_ids, obs_diffs, pred_diffs, instance_masks, n_bootstraps):
    """
    For each instance mask, among masked instances, compute the regression slope of observed
    differences on predicted differences (calibration slope; 1 = calibrated magnitudes).
    Standard errors come from a gene-level block bootstrap, propagated through per-gene sufficient
    statistics (n, sum_x, sum_y, sum_xx, sum_xy with x=predicted, y=observed).
    """
    unique_genes, gene_indices = np.unique(gene_ids, return_inverse=True)
    n_genes = len(unique_genes)
    n_masks = len(instance_masks)

    gene_ns = np.zeros((n_genes, n_masks))
    gene_sum_xs = np.zeros((n_genes, n_masks))
    gene_sum_ys = np.zeros((n_genes, n_masks))
    gene_sum_xxs = np.zeros((n_genes, n_masks))
    gene_sum_xys = np.zeros((n_genes, n_masks))
    for mask_iter, instance_mask in enumerate(instance_masks):
        masked_gene_indices = gene_indices[instance_mask]
        gene_ns[:, mask_iter] = np.bincount(masked_gene_indices, minlength=n_genes)
        gene_sum_xs[:, mask_iter] = np.bincount(masked_gene_indices, weights=pred_diffs[instance_mask], minlength=n_genes)
        gene_sum_ys[:, mask_iter] = np.bincount(masked_gene_indices, weights=obs_diffs[instance_mask], minlength=n_genes)
        gene_sum_xxs[:, mask_iter] = np.bincount(masked_gene_indices, weights=np.square(pred_diffs[instance_mask]), minlength=n_genes)
        gene_sum_xys[:, mask_iter] = np.bincount(masked_gene_indices, weights=pred_diffs[instance_mask]*obs_diffs[instance_mask], minlength=n_genes)

    def slope_from_sufficient_stats(nn, sum_x, sum_y, sum_xx, sum_xy):
        with np.errstate(invalid='ignore', divide='ignore'):
            return (sum_xy - (sum_x*sum_y/nn))/(sum_xx - (np.square(sum_x)/nn))

    obs_vs_pred_slopes = slope_from_sufficient_stats(np.sum(gene_ns, axis=0), np.sum(gene_sum_xs, axis=0), np.sum(gene_sum_ys, axis=0), np.sum(gene_sum_xxs, axis=0), np.sum(gene_sum_xys, axis=0))

    # Gene-level block bootstrap
    np.random.seed(0)
    bs_slopes = np.zeros((n_bootstraps, n_masks))
    for bs_iter in range(n_bootstraps):
        bs_gene_multiplicities = np.bincount(np.random.choice(n_genes, size=n_genes, replace=True), minlength=n_genes)
        bs_slopes[bs_iter, :] = slope_from_sufficient_stats(np.dot(bs_gene_multiplicities, gene_ns), np.dot(bs_gene_multiplicities, gene_sum_xs), np.dot(bs_gene_multiplicities, gene_sum_ys), np.dot(bs_gene_multiplicities, gene_sum_xxs), np.dot(bs_gene_multiplicities, gene_sum_xys))
    obs_vs_pred_slope_ses = np.nanstd(bs_slopes, axis=0)

    n_instances = np.sum(gene_ns, axis=0)

    return obs_vs_pred_slopes, obs_vs_pred_slope_ses, n_instances


def compute_obs_vs_pred_slope_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snrs, snr_thresholds, n_bootstraps):
    """
    At each threshold, among instances whose stratification value (the snrs argument; any
    quantity works) exceeds the threshold, compute the observed-vs-predicted regression slope.
    """
    instance_masks = [snrs > snr_threshold for snr_threshold in snr_thresholds]
    return compute_obs_vs_pred_slope_over_instance_masks(gene_ids, obs_diffs, pred_diffs, instance_masks, n_bootstraps)


def assign_quantile_bins(values, n_bins):
    """
    Split instances into n_bins equal-count (quantile) bins of values. Returns the per-instance
    bin assignment (0 = lowest bin) and the average value within each bin (the bin's x-coordinate
    for plotting).
    """
    bin_edges = np.quantile(values, np.linspace(0.0, 1.0, n_bins + 1))
    bin_assignments = np.searchsorted(bin_edges, values, side='right') - 1
    bin_assignments = np.clip(bin_assignments, 0, n_bins - 1)
    bin_sums = np.bincount(bin_assignments, weights=values, minlength=n_bins)
    bin_counts = np.bincount(bin_assignments, minlength=n_bins)
    with np.errstate(invalid='ignore', divide='ignore'):
        bin_avg_values = bin_sums/bin_counts
    return bin_assignments, bin_avg_values


def compute_over_quantile_bins_in_chunks(core_fn, gene_ids, obs_diffs, pred_diffs, bin_assignments, n_bins, n_bootstraps, bin_chunk_size=100):
    """
    Run a *_over_instance_masks core one chunk of bins at a time and concatenate the results.
    Caps peak memory: the cores allocate (n_genes x n_masks) matrices, so a full call at large
    n_bins is several GB while chunks stay small. Results match the unchunked call because each
    core call reseeds its bootstrap RNG, so every chunk sees the identical gene multiplicities.
    """
    chunk_results = []
    for chunk_start in range(0, n_bins, bin_chunk_size):
        chunk_bin_masks = [bin_assignments == bin_iter for bin_iter in range(chunk_start, min(chunk_start + bin_chunk_size, n_bins))]
        chunk_results.append(core_fn(gene_ids, obs_diffs, pred_diffs, chunk_bin_masks, n_bootstraps))
    return tuple(np.concatenate(result_arrays) for result_arrays in zip(*chunk_results))


def compute_binned_observed_vs_predicted_differences(gene_ids, obs_diffs, pred_diffs, n_pred_diff_bins, n_bootstraps):
    """
    Bin instances into signed quantile bins of the predicted difference. Per bin, compute the
    average predicted difference and the average observed difference (SEs on the observed averages
    from a gene-level block bootstrap). On-the-line behavior (avg observed == avg predicted per bin)
    indicates calibrated prediction magnitudes.
    """
    unique_genes, gene_indices = np.unique(gene_ids, return_inverse=True)
    n_genes = len(unique_genes)

    # Signed quantile bins of the predicted difference
    bin_edges = np.quantile(pred_diffs, np.linspace(0.0, 1.0, n_pred_diff_bins + 1))
    bin_assignments = np.searchsorted(bin_edges, pred_diffs, side='right') - 1
    bin_assignments = np.clip(bin_assignments, 0, n_pred_diff_bins - 1)

    # Per-gene sums and counts in each bin
    gene_bin_obs_sums = np.zeros((n_genes, n_pred_diff_bins))
    gene_bin_pred_sums = np.zeros((n_genes, n_pred_diff_bins))
    gene_bin_counts = np.zeros((n_genes, n_pred_diff_bins))
    for bin_iter in range(n_pred_diff_bins):
        bin_indices = bin_assignments == bin_iter
        gene_bin_obs_sums[:, bin_iter] = np.bincount(gene_indices[bin_indices], weights=obs_diffs[bin_indices], minlength=n_genes)
        gene_bin_pred_sums[:, bin_iter] = np.bincount(gene_indices[bin_indices], weights=pred_diffs[bin_indices], minlength=n_genes)
        gene_bin_counts[:, bin_iter] = np.bincount(gene_indices[bin_indices], minlength=n_genes)

    with np.errstate(invalid='ignore', divide='ignore'):
        avg_pred_diffs = np.sum(gene_bin_pred_sums, axis=0)/np.sum(gene_bin_counts, axis=0)
        avg_obs_diffs = np.sum(gene_bin_obs_sums, axis=0)/np.sum(gene_bin_counts, axis=0)

    # Gene-level block bootstrap (on the observed averages)
    np.random.seed(0)
    bs_avg_obs_diffs = np.zeros((n_bootstraps, n_pred_diff_bins))
    for bs_iter in range(n_bootstraps):
        bs_gene_multiplicities = np.bincount(np.random.choice(n_genes, size=n_genes, replace=True), minlength=n_genes)
        with np.errstate(invalid='ignore', divide='ignore'):
            bs_avg_obs_diffs[bs_iter, :] = np.dot(bs_gene_multiplicities, gene_bin_obs_sums)/np.dot(bs_gene_multiplicities, gene_bin_counts)
    avg_obs_diff_ses = np.nanstd(bs_avg_obs_diffs, axis=0)

    n_instances = np.sum(gene_bin_counts, axis=0)

    return avg_pred_diffs, avg_obs_diffs, avg_obs_diff_ses, n_instances


def compute_per_tissue_pair_observed_differences_at_snr_threshold(gene_ids, obs_diffs, pred_diffs, snrs, instance_tissue_pair_labels, tissue_pair_labels, snr_threshold, n_bootstraps):
    """
    For each tissue pair separately, the average observed difference among predicted
    positives/negatives at a single fixed SNR threshold (SEs from gene-level block bootstrap
    within the tissue pair).
    """
    per_pair_results = []
    for tissue_pair_label in tissue_pair_labels:
        pair_indices = instance_tissue_pair_labels == tissue_pair_label
        avg_pos, avg_pos_se, avg_neg, avg_neg_se, n_pos, n_neg = compute_observed_differences_stratified_by_snr_thresholds(gene_ids[pair_indices], obs_diffs[pair_indices], pred_diffs[pair_indices], snrs[pair_indices], np.asarray([snr_threshold]), n_bootstraps)
        per_pair_results.append((tissue_pair_label, avg_pos[0], avg_pos_se[0], avg_neg[0], avg_neg_se[0], n_pos[0], n_neg[0]))
    return per_pair_results


parser = argparse.ArgumentParser(description='Merge per-tissue-pair expression differences results.')
parser.add_argument('--expression-differences-results-dir', dest='expression_differences_results_dir', required=True, help='Directory containing the per-tissue-pair expression differences summary files.')
parser.add_argument('--anno-method', dest='anno_method', required=True, help='Name of the annotation method used to generate the borzoi annotation files.')
args = parser.parse_args()

expression_differences_results_dir = args.expression_differences_results_dir
anno_method = args.anno_method

tissue_sample_pairs=["Heart_Left_Ventricle:GTEX-18465-0926-SM-731AY.1", "Brain_Cortex:GTEX-1H3O1-1726-SM-9WYSR.1", "Liver:GTEX-11EQ9-0526-SM-5A5JZ.1", "Whole_Blood:GTEX-1LB8K-0005-SM-DIPED.1", "Muscle_Skeletal:GTEX-13QJ3-0726-SM-5SI68.1"]
n_bootstraps=500
n_pred_diff_bins=20
n_stratification_quantile_bins=10000
per_tissue_pair_snr_threshold=0.5


# File names that the difference analysis was run on (one per unordered tissue pair)
expression_difference_result_files, tissue_pair_labels = extract_expression_difference_result_files(expression_differences_results_dir, tissue_sample_pairs, anno_method)


# Load in results across all tissue pairs into a single data structure
expression_differnce_results = load_in_expression_differences_results(expression_difference_result_files, tissue_pair_labels)
gene_ids, obs_diffs, pred_diffs, snrs, snr_alt_taus, instance_tissue_pair_labels = organize_expression_difference_results_into_arrays(expression_differnce_results)


# Continuous grid of SNR thresholds (0 through the 99.99th percentile of observed SNRs)
snr_thresholds = np.linspace(0.0, np.quantile(snrs, 0.9999), 200)

# Continuous grid of alternative-tau SNR thresholds (nan-aware: instances with undefined
# alternative-tau SNR carry nan and are excluded by every threshold mask)
snr_alt_tau_thresholds = np.linspace(0.0, np.nanquantile(snr_alt_taus, 0.9999), 200)

# Continuous grid of absolute predicted difference thresholds (0 through the 99.99th percentile)
abs_pred_diffs = np.abs(pred_diffs)
abs_pred_diff_thresholds = np.linspace(0.0, np.quantile(abs_pred_diffs, 0.9999), 200)


##################################
# 1. Average observed difference among predicted positives/negatives at each SNR threshold (SEs from gene-level block bootstrap)
##################################
avg_obs_diff_pos, avg_obs_diff_pos_se, avg_obs_diff_neg, avg_obs_diff_neg_se, n_pos_instances, n_neg_instances = compute_observed_differences_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snrs, snr_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'observed_differences_stratified_by_snr_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('snr_threshold\tavg_observed_diff_predicted_positives\tavg_observed_diff_predicted_positives_se\tavg_observed_diff_predicted_negatives\tavg_observed_diff_predicted_negatives_se\tn_predicted_positives\tn_predicted_negatives\n')
for threshold_iter, snr_threshold in enumerate(snr_thresholds):
    t.write(str(snr_threshold) + '\t' + str(avg_obs_diff_pos[threshold_iter]) + '\t' + str(avg_obs_diff_pos_se[threshold_iter]))
    t.write('\t' + str(avg_obs_diff_neg[threshold_iter]) + '\t' + str(avg_obs_diff_neg_se[threshold_iter]))
    t.write('\t' + str(int(n_pos_instances[threshold_iter])) + '\t' + str(int(n_neg_instances[threshold_iter])) + '\n')
t.close()
print(output_file)

# Same analysis stratified by the alternative-tau SNR
at_avg_obs_diff_pos, at_avg_obs_diff_pos_se, at_avg_obs_diff_neg, at_avg_obs_diff_neg_se, at_n_pos_instances, at_n_neg_instances = compute_observed_differences_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snr_alt_taus, snr_alt_tau_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'observed_differences_stratified_by_snr_alt_tau_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('snr_threshold\tavg_observed_diff_predicted_positives\tavg_observed_diff_predicted_positives_se\tavg_observed_diff_predicted_negatives\tavg_observed_diff_predicted_negatives_se\tn_predicted_positives\tn_predicted_negatives\n')
for threshold_iter, snr_threshold in enumerate(snr_alt_tau_thresholds):
    t.write(str(snr_threshold) + '\t' + str(at_avg_obs_diff_pos[threshold_iter]) + '\t' + str(at_avg_obs_diff_pos_se[threshold_iter]))
    t.write('\t' + str(at_avg_obs_diff_neg[threshold_iter]) + '\t' + str(at_avg_obs_diff_neg_se[threshold_iter]))
    t.write('\t' + str(int(at_n_pos_instances[threshold_iter])) + '\t' + str(int(at_n_neg_instances[threshold_iter])) + '\n')
t.close()
print(output_file)


##################################
# 2. Sign concordance at each SNR threshold
##################################
sign_concordances, sign_concordance_ses, n_concordance_instances = compute_sign_concordance_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snrs, snr_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'sign_concordance_stratified_by_snr_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('snr_threshold\tsign_concordance\tsign_concordance_se\tn_instances\n')
for threshold_iter, snr_threshold in enumerate(snr_thresholds):
    t.write(str(snr_threshold) + '\t' + str(sign_concordances[threshold_iter]) + '\t' + str(sign_concordance_ses[threshold_iter]) + '\t' + str(int(n_concordance_instances[threshold_iter])) + '\n')
t.close()
print(output_file)

# Same analysis stratified by the alternative-tau SNR
at_sign_concordances, at_sign_concordance_ses, at_n_concordance_instances = compute_sign_concordance_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snr_alt_taus, snr_alt_tau_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'sign_concordance_stratified_by_snr_alt_tau_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('snr_threshold\tsign_concordance\tsign_concordance_se\tn_instances\n')
for threshold_iter, snr_threshold in enumerate(snr_alt_tau_thresholds):
    t.write(str(snr_threshold) + '\t' + str(at_sign_concordances[threshold_iter]) + '\t' + str(at_sign_concordance_ses[threshold_iter]) + '\t' + str(int(at_n_concordance_instances[threshold_iter])) + '\n')
t.close()
print(output_file)


##################################
# 3. Regression slope of observed on predicted differences at each SNR threshold
##################################
obs_vs_pred_slopes, obs_vs_pred_slope_ses, n_slope_instances = compute_obs_vs_pred_slope_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snrs, snr_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'observed_vs_predicted_slope_stratified_by_snr_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('snr_threshold\tobs_vs_pred_slope\tobs_vs_pred_slope_se\tn_instances\n')
for threshold_iter, snr_threshold in enumerate(snr_thresholds):
    t.write(str(snr_threshold) + '\t' + str(obs_vs_pred_slopes[threshold_iter]) + '\t' + str(obs_vs_pred_slope_ses[threshold_iter]) + '\t' + str(int(n_slope_instances[threshold_iter])) + '\n')
t.close()
print(output_file)

# Same analysis stratified by the alternative-tau SNR
at_obs_vs_pred_slopes, at_obs_vs_pred_slope_ses, at_n_slope_instances = compute_obs_vs_pred_slope_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, snr_alt_taus, snr_alt_tau_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'observed_vs_predicted_slope_stratified_by_snr_alt_tau_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('snr_threshold\tobs_vs_pred_slope\tobs_vs_pred_slope_se\tn_instances\n')
for threshold_iter, snr_threshold in enumerate(snr_alt_tau_thresholds):
    t.write(str(snr_threshold) + '\t' + str(at_obs_vs_pred_slopes[threshold_iter]) + '\t' + str(at_obs_vs_pred_slope_ses[threshold_iter]) + '\t' + str(int(at_n_slope_instances[threshold_iter])) + '\n')
t.close()
print(output_file)


##################################
# 4. Average observed vs average predicted difference in signed quantile bins of the predicted difference
##################################
avg_pred_diffs_binned, avg_obs_diffs_binned, avg_obs_diff_binned_ses, n_binned_instances = compute_binned_observed_vs_predicted_differences(gene_ids, obs_diffs, pred_diffs, n_pred_diff_bins, n_bootstraps)

output_file = expression_differences_results_dir + 'observed_vs_predicted_difference_calibration_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('bin_index\tavg_predicted_diff\tavg_observed_diff\tavg_observed_diff_se\tn_instances\n')
for bin_iter in range(n_pred_diff_bins):
    t.write(str(bin_iter) + '\t' + str(avg_pred_diffs_binned[bin_iter]) + '\t' + str(avg_obs_diffs_binned[bin_iter]) + '\t' + str(avg_obs_diff_binned_ses[bin_iter]) + '\t' + str(int(n_binned_instances[bin_iter])) + '\n')
t.close()
print(output_file)


##################################
# 5. Per-tissue-pair observed differences at a fixed SNR threshold
##################################
per_pair_results = compute_per_tissue_pair_observed_differences_at_snr_threshold(gene_ids, obs_diffs, pred_diffs, snrs, instance_tissue_pair_labels, tissue_pair_labels, per_tissue_pair_snr_threshold, n_bootstraps)

output_file = expression_differences_results_dir + 'per_tissue_pair_observed_differences_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('tissue_pair\tsnr_threshold\tavg_observed_diff_predicted_positives\tavg_observed_diff_predicted_positives_se\tavg_observed_diff_predicted_negatives\tavg_observed_diff_predicted_negatives_se\tn_predicted_positives\tn_predicted_negatives\n')
for per_pair_result in per_pair_results:
    t.write(per_pair_result[0] + '\t' + str(per_tissue_pair_snr_threshold) + '\t' + str(per_pair_result[1]) + '\t' + str(per_pair_result[2]))
    t.write('\t' + str(per_pair_result[3]) + '\t' + str(per_pair_result[4]) + '\t' + str(int(per_pair_result[5])) + '\t' + str(int(per_pair_result[6])) + '\n')
t.close()
print(output_file)

# Same analysis thresholding on the alternative-tau SNR
at_per_pair_results = compute_per_tissue_pair_observed_differences_at_snr_threshold(gene_ids, obs_diffs, pred_diffs, snr_alt_taus, instance_tissue_pair_labels, tissue_pair_labels, per_tissue_pair_snr_threshold, n_bootstraps)

output_file = expression_differences_results_dir + 'per_tissue_pair_observed_differences_snr_alt_tau_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('tissue_pair\tsnr_threshold\tavg_observed_diff_predicted_positives\tavg_observed_diff_predicted_positives_se\tavg_observed_diff_predicted_negatives\tavg_observed_diff_predicted_negatives_se\tn_predicted_positives\tn_predicted_negatives\n')
for per_pair_result in at_per_pair_results:
    t.write(per_pair_result[0] + '\t' + str(per_tissue_pair_snr_threshold) + '\t' + str(per_pair_result[1]) + '\t' + str(per_pair_result[2]))
    t.write('\t' + str(per_pair_result[3]) + '\t' + str(per_pair_result[4]) + '\t' + str(int(per_pair_result[5])) + '\t' + str(int(per_pair_result[6])) + '\n')
t.close()
print(output_file)


##################################
# 6. Average observed difference among predicted positives/negatives at each absolute predicted difference threshold
# (the SNR-stratification functions threshold on whatever values occupy the snr slot)
##################################
apd_avg_obs_diff_pos, apd_avg_obs_diff_pos_se, apd_avg_obs_diff_neg, apd_avg_obs_diff_neg_se, apd_n_pos_instances, apd_n_neg_instances = compute_observed_differences_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, abs_pred_diffs, abs_pred_diff_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'observed_differences_stratified_by_abs_pred_diff_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('abs_pred_diff_threshold\tavg_observed_diff_predicted_positives\tavg_observed_diff_predicted_positives_se\tavg_observed_diff_predicted_negatives\tavg_observed_diff_predicted_negatives_se\tn_predicted_positives\tn_predicted_negatives\n')
for threshold_iter, abs_pred_diff_threshold in enumerate(abs_pred_diff_thresholds):
    t.write(str(abs_pred_diff_threshold) + '\t' + str(apd_avg_obs_diff_pos[threshold_iter]) + '\t' + str(apd_avg_obs_diff_pos_se[threshold_iter]))
    t.write('\t' + str(apd_avg_obs_diff_neg[threshold_iter]) + '\t' + str(apd_avg_obs_diff_neg_se[threshold_iter]))
    t.write('\t' + str(int(apd_n_pos_instances[threshold_iter])) + '\t' + str(int(apd_n_neg_instances[threshold_iter])) + '\n')
t.close()
print(output_file)


##################################
# 7. Sign concordance at each absolute predicted difference threshold
##################################
apd_sign_concordances, apd_sign_concordance_ses, apd_n_concordance_instances = compute_sign_concordance_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, abs_pred_diffs, abs_pred_diff_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'sign_concordance_stratified_by_abs_pred_diff_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('abs_pred_diff_threshold\tsign_concordance\tsign_concordance_se\tn_instances\n')
for threshold_iter, abs_pred_diff_threshold in enumerate(abs_pred_diff_thresholds):
    t.write(str(abs_pred_diff_threshold) + '\t' + str(apd_sign_concordances[threshold_iter]) + '\t' + str(apd_sign_concordance_ses[threshold_iter]) + '\t' + str(int(apd_n_concordance_instances[threshold_iter])) + '\n')
t.close()
print(output_file)


##################################
# 8. Regression slope of observed on predicted differences at each absolute predicted difference threshold
##################################
apd_obs_vs_pred_slopes, apd_obs_vs_pred_slope_ses, apd_n_slope_instances = compute_obs_vs_pred_slope_stratified_by_snr_thresholds(gene_ids, obs_diffs, pred_diffs, abs_pred_diffs, abs_pred_diff_thresholds, n_bootstraps)

output_file = expression_differences_results_dir + 'observed_vs_predicted_slope_stratified_by_abs_pred_diff_thresholds_' + anno_method + '.txt'
t = open(output_file, 'w')
t.write('abs_pred_diff_threshold\tobs_vs_pred_slope\tobs_vs_pred_slope_se\tn_instances\n')
for threshold_iter, abs_pred_diff_threshold in enumerate(abs_pred_diff_thresholds):
    t.write(str(abs_pred_diff_threshold) + '\t' + str(apd_obs_vs_pred_slopes[threshold_iter]) + '\t' + str(apd_obs_vs_pred_slope_ses[threshold_iter]) + '\t' + str(int(apd_n_slope_instances[threshold_iter])) + '\n')
t.close()
print(output_file)


##################################
# 9. Observed differences, sign concordance, and obs-vs-pred slope in equal-count quantile bins
# of SNR and (separately) of the absolute predicted difference. Unlike the threshold analyses
# above (nested tails), bins are non-overlapping and each holds ~1/n_bins of the instances.
##################################
# Instances with an undefined alternative-tau SNR (nan) cannot be quantile-binned and are dropped
# from that stratification only
defined_alt_tau_indices = np.isfinite(snr_alt_taus)

for stratification_name, stratification_values, strat_gene_ids, strat_obs_diffs, strat_pred_diffs in [('snr', snrs, gene_ids, obs_diffs, pred_diffs), ('abs_pred_diff', abs_pred_diffs, gene_ids, obs_diffs, pred_diffs), ('snr_alt_tau', snr_alt_taus[defined_alt_tau_indices], gene_ids[defined_alt_tau_indices], obs_diffs[defined_alt_tau_indices], pred_diffs[defined_alt_tau_indices])]:
    bin_assignments, bin_avg_values = assign_quantile_bins(stratification_values, n_stratification_quantile_bins)

    bin_avg_obs_diff_pos, bin_avg_obs_diff_pos_se, bin_avg_obs_diff_neg, bin_avg_obs_diff_neg_se, bin_n_pos_instances, bin_n_neg_instances = compute_over_quantile_bins_in_chunks(compute_observed_differences_over_instance_masks, strat_gene_ids, strat_obs_diffs, strat_pred_diffs, bin_assignments, n_stratification_quantile_bins, n_bootstraps)
    output_file = expression_differences_results_dir + 'observed_differences_in_' + stratification_name + '_quantile_bins_' + anno_method + '.txt'
    t = open(output_file, 'w')
    t.write('bin_index\tavg_' + stratification_name + '\tavg_observed_diff_predicted_positives\tavg_observed_diff_predicted_positives_se\tavg_observed_diff_predicted_negatives\tavg_observed_diff_predicted_negatives_se\tn_predicted_positives\tn_predicted_negatives\n')
    for bin_iter in range(n_stratification_quantile_bins):
        t.write(str(bin_iter) + '\t' + str(bin_avg_values[bin_iter]) + '\t' + str(bin_avg_obs_diff_pos[bin_iter]) + '\t' + str(bin_avg_obs_diff_pos_se[bin_iter]))
        t.write('\t' + str(bin_avg_obs_diff_neg[bin_iter]) + '\t' + str(bin_avg_obs_diff_neg_se[bin_iter]))
        t.write('\t' + str(int(bin_n_pos_instances[bin_iter])) + '\t' + str(int(bin_n_neg_instances[bin_iter])) + '\n')
    t.close()
    print(output_file)

    bin_sign_concordances, bin_sign_concordance_ses, bin_n_concordance_instances = compute_over_quantile_bins_in_chunks(compute_sign_concordance_over_instance_masks, strat_gene_ids, strat_obs_diffs, strat_pred_diffs, bin_assignments, n_stratification_quantile_bins, n_bootstraps)
    output_file = expression_differences_results_dir + 'sign_concordance_in_' + stratification_name + '_quantile_bins_' + anno_method + '.txt'
    t = open(output_file, 'w')
    t.write('bin_index\tavg_' + stratification_name + '\tsign_concordance\tsign_concordance_se\tn_instances\n')
    for bin_iter in range(n_stratification_quantile_bins):
        t.write(str(bin_iter) + '\t' + str(bin_avg_values[bin_iter]) + '\t' + str(bin_sign_concordances[bin_iter]) + '\t' + str(bin_sign_concordance_ses[bin_iter]) + '\t' + str(int(bin_n_concordance_instances[bin_iter])) + '\n')
    t.close()
    print(output_file)

    bin_obs_vs_pred_slopes, bin_obs_vs_pred_slope_ses, bin_n_slope_instances = compute_over_quantile_bins_in_chunks(compute_obs_vs_pred_slope_over_instance_masks, strat_gene_ids, strat_obs_diffs, strat_pred_diffs, bin_assignments, n_stratification_quantile_bins, n_bootstraps)
    output_file = expression_differences_results_dir + 'observed_vs_predicted_slope_in_' + stratification_name + '_quantile_bins_' + anno_method + '.txt'
    t = open(output_file, 'w')
    t.write('bin_index\tavg_' + stratification_name + '\tobs_vs_pred_slope\tobs_vs_pred_slope_se\tn_instances\n')
    for bin_iter in range(n_stratification_quantile_bins):
        t.write(str(bin_iter) + '\t' + str(bin_avg_values[bin_iter]) + '\t' + str(bin_obs_vs_pred_slopes[bin_iter]) + '\t' + str(bin_obs_vs_pred_slope_ses[bin_iter]) + '\t' + str(int(bin_n_slope_instances[bin_iter])) + '\n')
    t.close()
    print(output_file)
