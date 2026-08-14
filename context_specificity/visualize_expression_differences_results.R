args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(scales)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_observed_difference_by_snr_threshold_plot <- function(df, plot_title) {
	# Line plot of the average observed expression difference among instances passing each SNR
	# threshold, separately for predicted positives (pred diff > 0) and predicted negatives
	# (pred diff < 0). Ribbons are gaussian-approximation 95% CIs from the gene-block bootstrap SEs.
	long_df = rbind(
		data.frame(
			snr_threshold=df$snr_threshold,
			prediction_class="Predicted positives",
			avg_observed_diff=df$avg_observed_diff_predicted_positives,
			se=df$avg_observed_diff_predicted_positives_se,
			stringsAsFactors=FALSE
		),
		data.frame(
			snr_threshold=df$snr_threshold,
			prediction_class="Predicted negatives",
			avg_observed_diff=df$avg_observed_diff_predicted_negatives,
			se=df$avg_observed_diff_predicted_negatives_se,
			stringsAsFactors=FALSE
		)
	)
	long_df = long_df[is.finite(long_df$avg_observed_diff), ]
	long_df$prediction_class = factor(long_df$prediction_class, levels=c("Predicted positives", "Predicted negatives"))
	# Gaussian approximation confidence intervals
	long_df$ci_lower = long_df$avg_observed_diff - 1.96*long_df$se
	long_df$ci_upper = long_df$avg_observed_diff + 1.96*long_df$se
	return(
		ggplot(long_df, aes(x=snr_threshold, y=avg_observed_diff, color=prediction_class, fill=prediction_class)) +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=.18, color=NA) +
		geom_line(linewidth=.7) +
		scale_color_manual(values=c("Predicted positives"="#3E7A34", "Predicted negatives"="#D06A4B"), name="") +
		scale_fill_manual(values=c("Predicted positives"="#3E7A34", "Predicted negatives"="#D06A4B"), name="") +
		xlab("SNR threshold") +
		ylab("Average observed expression difference") +
		ggtitle(plot_title) +
		figure_theme() +
		theme(legend.position="right")
	)
}


make_observed_vs_predicted_calibration_plot <- function(df, plot_title) {
	# Average observed vs average predicted difference in signed quantile bins of the predicted
	# difference. On-the-line (y=x) behavior indicates calibrated prediction magnitudes.
	df = df[is.finite(df$avg_observed_diff), ]
	# Gaussian approximation confidence intervals
	df$ci_lower = df$avg_observed_diff - 1.96*df$avg_observed_diff_se
	df$ci_upper = df$avg_observed_diff + 1.96*df$avg_observed_diff_se
	return(
		ggplot(df, aes(x=avg_predicted_diff, y=avg_observed_diff)) +
		geom_abline(intercept=0, slope=1, linetype="dashed", linewidth=.4, color="#6B7280") +
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=0, linewidth=.45, color="#33424F") +
		geom_point(size=2, color="#4C78A8") +
		xlab("Average predicted expression difference") +
		ylab("Average observed expression difference") +
		ggtitle(plot_title) +
		figure_theme()
	)
}


make_sign_concordance_plot <- function(df, plot_title) {
	# Fraction of instances passing each SNR threshold where sign(observed diff) == sign(predicted diff).
	df = df[is.finite(df$sign_concordance), ]
	# Gaussian approximation confidence intervals
	df$ci_lower = df$sign_concordance - 1.96*df$sign_concordance_se
	df$ci_upper = df$sign_concordance + 1.96*df$sign_concordance_se
	return(
		ggplot(df, aes(x=snr_threshold, y=sign_concordance)) +
		geom_hline(yintercept=0.5, linewidth=.4, color="#6B7280", linetype="dashed") +
		geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=.18, fill="#4C78A8") +
		geom_line(linewidth=.7, color="#4C78A8") +
		xlab("SNR threshold") +
		ylab("Sign concordance") +
		ggtitle(plot_title) +
		figure_theme()
	)
}


make_obs_vs_pred_slope_plot <- function(df, plot_title) {
	# Regression slope of observed on predicted differences among instances passing each SNR
	# threshold. 1 (dashed line) indicates calibrated prediction magnitudes.
	df = df[is.finite(df$obs_vs_pred_slope), ]
	# Gaussian approximation confidence intervals
	df$ci_lower = df$obs_vs_pred_slope - 1.96*df$obs_vs_pred_slope_se
	df$ci_upper = df$obs_vs_pred_slope + 1.96*df$obs_vs_pred_slope_se
	return(
		ggplot(df, aes(x=snr_threshold, y=obs_vs_pred_slope)) +
		geom_hline(yintercept=1.0, linewidth=.4, color="#D06A4B", linetype="dashed") +
		geom_hline(yintercept=0.0, linewidth=.4, color="#6B7280", linetype="dashed") +
		geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=.18, fill="#4C78A8") +
		geom_line(linewidth=.7, color="#4C78A8") +
		xlab("SNR threshold") +
		ylab("Observed vs predicted slope") +
		ggtitle(plot_title) +
		figure_theme()
	)
}


make_instance_counts_plot <- function(df, plot_title) {
	# Number of instances passing each SNR threshold, per prediction class (log10 y-axis).
	long_df = rbind(
		data.frame(snr_threshold=df$snr_threshold, prediction_class="Predicted positives", n_instances=df$n_predicted_positives, stringsAsFactors=FALSE),
		data.frame(snr_threshold=df$snr_threshold, prediction_class="Predicted negatives", n_instances=df$n_predicted_negatives, stringsAsFactors=FALSE)
	)
	long_df = long_df[long_df$n_instances > 0, ]
	long_df$prediction_class = factor(long_df$prediction_class, levels=c("Predicted positives", "Predicted negatives"))
	return(
		ggplot(long_df, aes(x=snr_threshold, y=n_instances, color=prediction_class)) +
		geom_line(linewidth=.7) +
		scale_color_manual(values=c("Predicted positives"="#3E7A34", "Predicted negatives"="#D06A4B"), name="") +
		scale_y_log10(labels=comma_format()) +
		xlab("SNR threshold") +
		ylab("Instances passing threshold") +
		ggtitle(plot_title) +
		figure_theme() +
		theme(legend.position="top")
	)
}


make_tissue_pair_forest_plot <- function(df, plot_title) {
	# Per-tissue-pair average observed difference among predicted positives/negatives at the fixed
	# SNR threshold recorded in the file.
	long_df = rbind(
		data.frame(
			tissue_pair=df$tissue_pair,
			prediction_class="Predicted positives",
			avg_observed_diff=df$avg_observed_diff_predicted_positives,
			se=df$avg_observed_diff_predicted_positives_se,
			stringsAsFactors=FALSE
		),
		data.frame(
			tissue_pair=df$tissue_pair,
			prediction_class="Predicted negatives",
			avg_observed_diff=df$avg_observed_diff_predicted_negatives,
			se=df$avg_observed_diff_predicted_negatives_se,
			stringsAsFactors=FALSE
		)
	)
	long_df = long_df[is.finite(long_df$avg_observed_diff), ]
	long_df$prediction_class = factor(long_df$prediction_class, levels=c("Predicted positives", "Predicted negatives"))
	long_df$tissue_pair_label = gsub("_", " ", long_df$tissue_pair)
	long_df$tissue_pair_label = factor(long_df$tissue_pair_label, levels=rev(unique(long_df$tissue_pair_label)))
	# Gaussian approximation confidence intervals
	long_df$ci_lower = long_df$avg_observed_diff - 1.96*long_df$se
	long_df$ci_upper = long_df$avg_observed_diff + 1.96*long_df$se
	dodge_width <- 0.55
	return(
		ggplot(long_df, aes(x=tissue_pair_label, y=avg_observed_diff, color=prediction_class)) +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.18, linewidth=.45, position=position_dodge(width=dodge_width)) +
		geom_point(size=2, position=position_dodge(width=dodge_width)) +
		scale_color_manual(values=c("Predicted positives"="#3E7A34", "Predicted negatives"="#D06A4B"), name="") +
		coord_flip() +
		xlab("") +
		ylab("Average observed expression difference") +
		ggtitle(plot_title) +
		figure_theme() +
		theme(legend.position="top")
	)
}


#########################
# Command line args
#########################
expression_differences_results_dir = args[1]
anno_method = args[2]
visualization_dir = args[3]


#########################
# Load in merged expression differences results
#########################
merged_results_file = paste0(expression_differences_results_dir, "observed_differences_stratified_by_snr_thresholds_", anno_method, ".txt")
merged_df = read.table(merged_results_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)


#########################
# Observed expression differences stratified by SNR threshold
#########################
snr_threshold_plot = make_observed_difference_by_snr_threshold_plot(merged_df, "")
snr_threshold_plot_output_file = paste0(visualization_dir, "observed_differences_stratified_by_snr_thresholds_", anno_method, ".pdf")
ggsave(snr_threshold_plot_output_file, snr_threshold_plot + theme(legend.position="top"), width=7.2, height=3.6)


#########################
# Instances passing each SNR threshold (from the same merged file)
#########################
instance_counts_plot = make_instance_counts_plot(merged_df, "")
instance_counts_plot_output_file = paste0(visualization_dir, "instance_counts_stratified_by_snr_thresholds_", anno_method, ".pdf")
ggsave(instance_counts_plot_output_file, instance_counts_plot, width=7.2, height=3.6)


#########################
# Sign concordance stratified by SNR threshold
#########################
concordance_file = paste0(expression_differences_results_dir, "sign_concordance_stratified_by_snr_thresholds_", anno_method, ".txt")
concordance_df = read.table(concordance_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
sign_concordance_plot = make_sign_concordance_plot(concordance_df, "")
sign_concordance_plot_output_file = paste0(visualization_dir, "sign_concordance_stratified_by_snr_thresholds_", anno_method, ".pdf")
ggsave(sign_concordance_plot_output_file, sign_concordance_plot, width=7.2, height=3.6)


#########################
# Observed vs predicted slope stratified by SNR threshold
#########################
slope_file = paste0(expression_differences_results_dir, "observed_vs_predicted_slope_stratified_by_snr_thresholds_", anno_method, ".txt")
slope_df = read.table(slope_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
obs_vs_pred_slope_plot = make_obs_vs_pred_slope_plot(slope_df, "")
obs_vs_pred_slope_plot_output_file = paste0(visualization_dir, "observed_vs_predicted_slope_stratified_by_snr_thresholds_", anno_method, ".pdf")
ggsave(obs_vs_pred_slope_plot_output_file, obs_vs_pred_slope_plot, width=7.2, height=3.6)


#########################
# Observed vs predicted difference calibration (signed quantile bins of predicted difference)
#########################
calibration_file = paste0(expression_differences_results_dir, "observed_vs_predicted_difference_calibration_", anno_method, ".txt")
calibration_df = read.table(calibration_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
calibration_plot = make_observed_vs_predicted_calibration_plot(calibration_df, "")
calibration_plot_output_file = paste0(visualization_dir, "observed_vs_predicted_difference_calibration_", anno_method, ".pdf")
ggsave(calibration_plot_output_file, calibration_plot, width=4.8, height=4.2)


#########################
# Per-tissue-pair observed differences at a fixed SNR threshold
#########################
per_pair_file = paste0(expression_differences_results_dir, "per_tissue_pair_observed_differences_", anno_method, ".txt")
per_pair_df = read.table(per_pair_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
per_pair_snr_threshold = per_pair_df$snr_threshold[1]
tissue_pair_forest_plot = make_tissue_pair_forest_plot(per_pair_df, paste0("SNR threshold = ", per_pair_snr_threshold))
tissue_pair_forest_plot_output_file = paste0(visualization_dir, "per_tissue_pair_observed_differences_", anno_method, ".pdf")
ggsave(tissue_pair_forest_plot_output_file, tissue_pair_forest_plot, width=7.2, height=4.2)
