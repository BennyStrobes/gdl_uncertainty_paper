args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
library(scales)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



extract_intercept_correlation <- function(stats_file, target_annotation_name="intercept", target_category_name="intercept") {
    # Pull one correlation (mean, bootstrap SE) out of an S-LDMC bootstrap_stats file: the intercept's
    # by default, or any (annotation, category) cell's (e.g. borzoi_magnitude_bins / magnitude_bin4)
    stats_df = read.table(stats_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    target_row = stats_df[stats_df$annotation_name == target_annotation_name & stats_df$category_name == target_category_name & stats_df$output_name == "correlation", ]
    if (nrow(target_row) != 1) {
        stop(paste("Expected exactly one", target_annotation_name, "/", target_category_name, "correlation row in", stats_file, "but found", nrow(target_row)))
    }
    return(c(target_row$mean, target_row$bootstrap_se))
}


load_in_tissue_permuted_results <- function(tissue_permuted_pairs_file, tissue_permuted_sldmc_results_output_dir, sldmc_results_output_dir, target_annotation_name="intercept", target_category_name="intercept") {
    # One row per tissue, holding one cell's correlation (with 95% Gaussian CI) from the standard
    # S-LDMC run (borzoi effects from the matched tissue, in sldmc_results_output_dir) and from the
    # tissue-permuted run (borzoi effects from the randomly drawn partner tissue, in
    # tissue_permuted_sldmc_results_output_dir). By default the intercept's correlation is pulled;
    # target_annotation_name / target_category_name select any other cell.
    tissue_permuted_pairs_df = read.table(tissue_permuted_pairs_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

    merged_df = data.frame()
    for (i in 1:nrow(tissue_permuted_pairs_df)) {
        tissue = tissue_permuted_pairs_df$GTEx_tissue[i]
        target_identifier = tissue_permuted_pairs_df$target_identifier[i]
        permuted_tissue = tissue_permuted_pairs_df$permuted_GTEx_tissue[i]
        permuted_target_identifier = tissue_permuted_pairs_df$permuted_target_identifier[i]

        standard_stats_file = paste0(sldmc_results_output_dir, "sldmc_results_", tissue, "_", target_identifier, "_default_bootstrap_stats.txt")
        permuted_stats_file = paste0(tissue_permuted_sldmc_results_output_dir, "sldmc_results_eqtl_", tissue, "_borzoi_", permuted_tissue, "_", permuted_target_identifier, "_default_bootstrap_stats.txt")
        if (file.exists(standard_stats_file) == FALSE) {
            print(paste("WARNING: standard S-LDMC stats file not found; skipping tissue:", standard_stats_file))
            next
        }
        if (file.exists(permuted_stats_file) == FALSE) {
            print(paste("WARNING: tissue-permuted S-LDMC stats file not found; skipping tissue:", permuted_stats_file))
            next
        }

        standard_stats = extract_intercept_correlation(standard_stats_file, target_annotation_name, target_category_name)
        permuted_stats = extract_intercept_correlation(permuted_stats_file, target_annotation_name, target_category_name)

        merged_df = rbind(merged_df, data.frame(
            GTEx_tissue=tissue,
            target_identifier=target_identifier,
            permuted_GTEx_tissue=permuted_tissue,
            permuted_target_identifier=permuted_target_identifier,
            standard_correlation=standard_stats[1],
            standard_correlation_se=standard_stats[2],
            standard_correlation_ci_lower=standard_stats[1] - 1.96*standard_stats[2],
            standard_correlation_ci_upper=standard_stats[1] + 1.96*standard_stats[2],
            permuted_correlation=permuted_stats[1],
            permuted_correlation_se=permuted_stats[2],
            permuted_correlation_ci_lower=permuted_stats[1] - 1.96*permuted_stats[2],
            permuted_correlation_ci_upper=permuted_stats[1] + 1.96*permuted_stats[2]
        ))
    }
    if (nrow(merged_df) == 0) {
        stop("No tissues with both standard and tissue-permuted S-LDMC results loaded")
    }
    return(merged_df)
}


load_in_tissue_permuted_magnitude_bin_results <- function(tissue_permuted_pairs_file, tissue_permuted_sldmc_results_output_dir, sldmc_results_output_dir) {
    # Long version of load_in_tissue_permuted_results across the five coarse borzoi magnitude bins:
    # one row per (tissue, magnitude bin) with the matched and permuted correlations.
    bin_df_list = list()
    for (magnitude_bin_id in 0:4) {
        bin_df = load_in_tissue_permuted_results(tissue_permuted_pairs_file, tissue_permuted_sldmc_results_output_dir, sldmc_results_output_dir, target_annotation_name="borzoi_magnitude_bins", target_category_name=paste0("magnitude_bin", magnitude_bin_id))
        bin_df$magnitude_bin_id = magnitude_bin_id
        bin_df_list[[magnitude_bin_id + 1]] = bin_df
    }
    return(do.call(rbind, bin_df_list))
}


make_matched_minus_permuted_difference_dot_plot <- function(df) {
    # Jittered dot plot of the per-tissue paired difference in correlation (matched - permuted), one
    # dot per tissue, across the coarse borzoi magnitude bins. Above each bin: how many tissues fall
    # above zero and a two-sided sign test (binomial) p-value.
    magnitude_bin_labels = c("<0.001", "0.001-0.01", "0.01-0.075", "0.075-0.2", ">=0.2")
    plot_df = df
    plot_df$difference = plot_df$standard_correlation - plot_df$permuted_correlation
    plot_df$magnitude_bin = factor(plot_df$magnitude_bin_id, levels=0:4, labels=magnitude_bin_labels)

    annotation_df = data.frame()
    for (bin_label in levels(plot_df$magnitude_bin)) {
        bin_differences = plot_df$difference[plot_df$magnitude_bin == bin_label & plot_df$difference != 0]
        n_positive = sum(bin_differences > 0)
        sign_test_pvalue = binom.test(n_positive, length(bin_differences))$p.value
        annotation_df = rbind(annotation_df, data.frame(
            magnitude_bin=factor(bin_label, levels=levels(plot_df$magnitude_bin)),
            label=paste0(n_positive, "/", length(bin_differences), " > 0\np=", signif(sign_test_pvalue, 2))
        ))
    }

    return(
        ggplot(plot_df, aes(x=magnitude_bin, y=difference)) +
        geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
        geom_jitter(position=position_jitter(width=.15, height=0, seed=1), shape=21, size=1.8, stroke=.35, fill="#3F88C5", color="#111827", alpha=.8) +
        geom_text(data=annotation_df, aes(x=magnitude_bin, y=Inf, label=label), vjust=1.2, size=2.9, color="#111827") +
        # Extra headroom keeps the sign-test labels clear of the points
        scale_y_continuous(labels=number_format(accuracy=.01), expand=expansion(mult=c(.05, .25))) +
        xlab("Borzoi magnitude bin") +
        ylab("S-LDMC correlation difference\n(matched - permuted tissue)") +
        figure_theme() +
        theme(
            axis.text.x=element_text(angle=35, hjust=1, vjust=1),
            plot.margin=margin(8, 14, 8, 8)
        )
    )
}


make_standard_vs_permuted_correlation_scatter <- function(df, xlab_title="S-LDMC correlation\n(permuted tissue)", ylab_title="S-LDMC correlation\n(matched tissue)") {
    # Scatter of a correlation from the standard S-LDMC run (y) against the
    # tissue-permuted run (x), one point per tissue with 95% Gaussian CI error bars in both directions
    # and a dashed y=x reference line. Axes share the same limits so y=x runs corner to corner.
    axis_limits = c(
        min(c(df$standard_correlation_ci_lower, df$permuted_correlation_ci_lower)),
        max(c(df$standard_correlation_ci_upper, df$permuted_correlation_ci_upper))
    )
    return(
        ggplot(df, aes(x=permuted_correlation, y=standard_correlation)) +
        geom_abline(intercept=0, slope=1, linewidth=.4, color="#6B7280", linetype="dashed") +
        geom_errorbar(aes(ymin=standard_correlation_ci_lower, ymax=standard_correlation_ci_upper), width=0, linewidth=.35, color="#9CA3AF") +
        geom_errorbarh(aes(xmin=permuted_correlation_ci_lower, xmax=permuted_correlation_ci_upper), height=0, linewidth=.35, color="#9CA3AF") +
        geom_point(shape=21, size=1.8, stroke=.35, fill="#3F88C5", color="#111827", alpha=.8) +
        scale_x_continuous(labels=number_format(accuracy=.01), limits=axis_limits) +
        scale_y_continuous(labels=number_format(accuracy=.01), limits=axis_limits) +
        coord_fixed() +
        xlab(xlab_title) +
        ylab(ylab_title) +
        figure_theme() +
        theme(plot.margin=margin(8, 14, 8, 8))
    )
}


#######################
# Command line arguments
#######################
if (length(args) != 4) {
    stop("Usage: Rscript visualize_tissue_permuted_sldmc_results.R <tissue_permuted_sldmc_results_output_dir> <tissue_permuted_pairs_file> <visualize_tissue_permuted_sldmc_results_dir> <sldmc_results_output_dir>")
}
tissue_permuted_sldmc_results_output_dir = args[1]
tissue_permuted_pairs_file = args[2]
visualize_tissue_permuted_sldmc_results_dir = args[3]
# Directory holding the standard (non-permuted) S-LDMC results
sldmc_results_output_dir = args[4]



# Load in tissue permuted results
tissue_permuted_sldmc_results_df = load_in_tissue_permuted_results(tissue_permuted_pairs_file, tissue_permuted_sldmc_results_output_dir, sldmc_results_output_dir)


########################
# Scatter of the intercept correlation: tissue-permuted run against the standard run (one point per
# tissue, 95% CIs in both directions, dashed y=x line)
########################
standard_vs_permuted_scatter_plot = make_standard_vs_permuted_correlation_scatter(tissue_permuted_sldmc_results_df)
standard_vs_permuted_scatter_output_file = paste0(visualize_tissue_permuted_sldmc_results_dir, "sldmc_standard_vs_tissue_permuted_intercept_correlation_scatter.pdf")
ggsave(standard_vs_permuted_scatter_output_file, standard_vs_permuted_scatter_plot, width=4.6, height=4.4)


########################
# Same scatter, but for the correlation of large-magnitude borzoi variants: the top coarse borzoi
# magnitude bin (borzoi_magnitude_bins / magnitude_bin4, i.e. magnitude >= 0.2)
########################
tissue_permuted_large_magnitude_df = load_in_tissue_permuted_results(tissue_permuted_pairs_file, tissue_permuted_sldmc_results_output_dir, sldmc_results_output_dir, target_annotation_name="borzoi_magnitude_bins", target_category_name="magnitude_bin4")
large_magnitude_scatter_plot = make_standard_vs_permuted_correlation_scatter(tissue_permuted_large_magnitude_df, xlab_title="S-LDMC correlation, magnitude >= 0.2\n(permuted tissue)", ylab_title="S-LDMC correlation, magnitude >= 0.2\n(matched tissue)")
large_magnitude_scatter_output_file = paste0(visualize_tissue_permuted_sldmc_results_dir, "sldmc_standard_vs_tissue_permuted_large_magnitude_correlation_scatter.pdf")
ggsave(large_magnitude_scatter_output_file, large_magnitude_scatter_plot, width=4.6, height=4.4)


########################
# Per-tissue paired difference in correlation (matched - permuted) across the coarse borzoi magnitude
# bins (one dot per tissue per bin, with a per-bin two-sided sign test)
########################
tissue_permuted_magnitude_bin_df = load_in_tissue_permuted_magnitude_bin_results(tissue_permuted_pairs_file, tissue_permuted_sldmc_results_output_dir, sldmc_results_output_dir)
matched_minus_permuted_difference_plot = make_matched_minus_permuted_difference_dot_plot(tissue_permuted_magnitude_bin_df)
matched_minus_permuted_difference_output_file = paste0(visualize_tissue_permuted_sldmc_results_dir, "sldmc_matched_minus_permuted_correlation_by_magnitude_bin.pdf")
ggsave(matched_minus_permuted_difference_output_file, matched_minus_permuted_difference_plot, width=7.2, height=3.6)


