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
