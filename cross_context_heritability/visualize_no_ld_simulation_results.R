args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(scales)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}





compute_true_decile_summary_df <- function(df, true_col, est_col, n_bins) {
	# Split genes into n_bins equal-count bins (deciles when n_bins=10) by the simulated (true) value in
	# true_col. Each row gives the bin's mean true value and the mean estimate (est_col) with a 95% CI
	# (mean +/- 1.96 * empirical SE across genes in the bin).
	df = df[!is.na(df[[true_col]]) & !is.na(df[[est_col]]), ]
	bin_id = ceiling(rank(df[[true_col]], ties.method="first") / nrow(df) * n_bins)

	bin_arr = c()
	n_gene_arr = c()
	mean_true_arr = c()
	mean_est_arr = c()
	mean_est_lb_arr = c()
	mean_est_ub_arr = c()
	for (bin_iter in 1:n_bins) {
		bin_df = df[bin_id == bin_iter, ]
		mean_est = mean(bin_df[[est_col]])
		se_est = sd(bin_df[[est_col]]) / sqrt(nrow(bin_df))

		bin_arr = c(bin_arr, bin_iter)
		n_gene_arr = c(n_gene_arr, nrow(bin_df))
		mean_true_arr = c(mean_true_arr, mean(bin_df[[true_col]]))
		mean_est_arr = c(mean_est_arr, mean_est)
		mean_est_lb_arr = c(mean_est_lb_arr, mean_est - 1.96*se_est)
		mean_est_ub_arr = c(mean_est_ub_arr, mean_est + 1.96*se_est)
	}

	summary_df = data.frame(bin=bin_arr, n_genes=n_gene_arr, mean_true=mean_true_arr, mean_est=mean_est_arr, mean_est_lb=mean_est_lb_arr, mean_est_ub=mean_est_ub_arr)
	return(summary_df)
}


make_true_vs_est_decile_plot <- function(summary_df, xlab, ylab) {
	# Calibration plot: mean estimate (points with 95% CI, joined by a line) against the mean true value in
	# each decile of the true value. The dashed line is y = x.
	reference_color = "#3F3F46"
	point_color = "#2563EB"

	pp = ggplot(summary_df, aes(x=mean_true, y=mean_est)) +
		geom_abline(slope=1, intercept=0, color=reference_color, linetype="dashed", linewidth=0.6) +
		geom_line(color=point_color, linewidth=0.7) +
		geom_pointrange(aes(ymin=mean_est_lb, ymax=mean_est_ub), color=point_color, linewidth=0.7) +
		geom_point(color=point_color, size=2.4) +
		figure_theme() +
		labs(x=xlab, y=ylab)
	return(pp)
}


compute_ci_coverage_summary_df <- function(df, true_col, est_col, se_col, coverage_levels) {
	# For each nominal level L, the CI is est +/- qnorm((1+L)/2) * se. Each row gives the fraction of genes
	# whose CI contains the true value (observed coverage), with a 95% binomial (normal approximation) CI.
	df = df[!is.na(df[[true_col]]) & !is.na(df[[est_col]]) & !is.na(df[[se_col]]), ]
	n_genes = nrow(df)

	level_arr = c()
	coverage_arr = c()
	coverage_lb_arr = c()
	coverage_ub_arr = c()
	for (level in coverage_levels) {
		z = qnorm((1.0 + level) / 2.0)
		ci_lb = df[[est_col]] - z*df[[se_col]]
		ci_ub = df[[est_col]] + z*df[[se_col]]
		covered = (df[[true_col]] >= ci_lb) & (df[[true_col]] <= ci_ub)
		coverage = mean(covered)
		coverage_se = sqrt(coverage*(1.0 - coverage) / n_genes)

		level_arr = c(level_arr, level)
		coverage_arr = c(coverage_arr, coverage)
		coverage_lb_arr = c(coverage_lb_arr, coverage - 1.96*coverage_se)
		coverage_ub_arr = c(coverage_ub_arr, coverage + 1.96*coverage_se)
	}

	summary_df = data.frame(level=level_arr, coverage=coverage_arr, coverage_lb=coverage_lb_arr, coverage_ub=coverage_ub_arr, n_genes=n_genes)
	return(summary_df)
}


make_ci_coverage_plot <- function(summary_df, xlab, ylab) {
	# Coverage plot: observed coverage (points with 95% CI, joined by a line) against the nominal CI level.
	# The dashed line is y = x.
	reference_color = "#3F3F46"
	point_color = "#2563EB"

	pp = ggplot(summary_df, aes(x=level, y=coverage)) +
		geom_abline(slope=1, intercept=0, color=reference_color, linetype="dashed", linewidth=0.6) +
		geom_line(color=point_color, linewidth=0.7) +
		geom_pointrange(aes(ymin=coverage_lb, ymax=coverage_ub), color=point_color, linewidth=0.7) +
		geom_point(color=point_color, size=2.4) +
		scale_x_continuous(breaks=summary_df$level, labels=percent(summary_df$level)) +
		scale_y_continuous(labels=percent) +
		figure_theme() +
		labs(x=xlab, y=ylab)
	return(pp)
}


make_pvalue_histogram <- function(df, pvalue_col, bin_width, xlab, ylab) {
	# Histogram of p-values in bins of width bin_width on [0,1]. The dashed line is the expected count per
	# bin under uniformly distributed p-values.
	reference_color = "#3F3F46"
	bar_color = "#2563EB"
	df = df[!is.na(df[[pvalue_col]]), ]
	expected_count = nrow(df) * bin_width

	pp = ggplot(df, aes(x=.data[[pvalue_col]])) +
		geom_histogram(breaks=seq(0, 1, by=bin_width), fill=bar_color, color="white", linewidth=0.3) +
		geom_hline(yintercept=expected_count, color=reference_color, linetype="dashed", linewidth=0.6) +
		scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
		figure_theme() +
		labs(x=xlab, y=ylab)
	return(pp)
}


compute_power_summary_df <- function(df, pvalue_col, pvalue_thresholds) {
	# For each p-value threshold, the fraction of genes with p-value below the threshold (power), with a
	# 95% binomial (normal approximation) CI.
	df = df[!is.na(df[[pvalue_col]]), ]
	n_genes = nrow(df)

	threshold_arr = c()
	power_arr = c()
	power_lb_arr = c()
	power_ub_arr = c()
	for (threshold in pvalue_thresholds) {
		power = mean(df[[pvalue_col]] < threshold)
		power_se = sqrt(power*(1.0 - power) / n_genes)

		threshold_arr = c(threshold_arr, threshold)
		power_arr = c(power_arr, power)
		power_lb_arr = c(power_lb_arr, power - 1.96*power_se)
		power_ub_arr = c(power_ub_arr, power + 1.96*power_se)
	}

	summary_df = data.frame(threshold=threshold_arr, power=power_arr, power_lb=power_lb_arr, power_ub=power_ub_arr, n_genes=n_genes)
	summary_df$threshold_label = factor(paste0("p < ", summary_df$threshold), levels=paste0("p < ", summary_df$threshold))
	return(summary_df)
}


make_power_bar_plot <- function(summary_df, xlab, ylab) {
	# Bar plot of power (with 95% CI) at each p-value threshold. The dashed line over each bar marks the
	# threshold itself (the expected fraction of p-values below it under the null).
	reference_color = "#3F3F46"
	bar_color = "#2563EB"
	bar_width = 0.6

	pp = ggplot(summary_df, aes(x=threshold_label, y=power)) +
		geom_col(fill=bar_color, width=bar_width) +
		geom_errorbar(aes(ymin=power_lb, ymax=power_ub), width=0.18, linewidth=0.35, color="#111827") +
		geom_segment(aes(x=as.numeric(threshold_label) - bar_width/2, xend=as.numeric(threshold_label) + bar_width/2, y=threshold, yend=threshold), color=reference_color, linetype="dashed", linewidth=0.6) +
		scale_y_continuous(labels=percent, expand=expansion(mult=c(0, 0.05))) +
		expand_limits(y=0) +
		figure_theme() +
		labs(x=xlab, y=ylab)
	return(pp)
}



#####################
# Command line args
#####################
simulation_dir = args[1]
sample_size = args[2]
simulation_visualization_dir = args[3]


print(simulation_visualization_dir)

#####################
# Load in simulation results
#####################
sim_results_file = paste0(simulation_dir, "no_ld_simulation_results_", sample_size, ".txt")
sim_df = read.table(sim_results_file, header=TRUE, sep="\t")

null_sim_results_file = paste0(simulation_dir, "no_ld_null_simulation_results_", sample_size, ".txt")
null_sim_df = read.table(null_sim_results_file, header=TRUE, sep="\t")


#####################
# Figure 1: Estimated vs simulated delta beta variance (binned by decile of simulated delta beta variance)
#####################
decile_summary_df = compute_true_decile_summary_df(sim_df, "true_delta_beta_var", "est_delta_beta_var", 10)
decile_plot = make_true_vs_est_decile_plot(decile_summary_df, "Mean simulated delta beta variance", "Mean estimated delta beta variance\n")
output_file = paste0(simulation_visualization_dir, "no_ld_simulation_", sample_size, "_est_vs_true_delta_beta_var_deciles.pdf")
ggsave(output_file, decile_plot, width=7.2, height=4.5, units="in")



#####################
# Figure 2: Coverage of confidence intervals for delta beta variance (non-null simulations)
#####################
coverage_levels = c(0.75, 0.80, 0.85, 0.90, 0.95)
coverage_summary_df = compute_ci_coverage_summary_df(sim_df, "true_delta_beta_var", "est_delta_beta_var", "est_delta_beta_var_se", coverage_levels)
coverage_plot = make_ci_coverage_plot(coverage_summary_df, "Nominal confidence interval level", "Observed coverage")
output_file = paste0(simulation_visualization_dir, "no_ld_simulation_", sample_size, "_delta_beta_var_ci_coverage.pdf")
ggsave(output_file, coverage_plot, width=7.2, height=4.5, units="in")




#####################
# Figure 3: Histogram of p-values (null simulations)
#####################
null_pvalue_hist = make_pvalue_histogram(null_sim_df, "est_delta_beta_var_p", 0.05, "p-value", "Number of simulated genes")
output_file = paste0(simulation_visualization_dir, "no_ld_null_simulation_", sample_size, "_delta_beta_var_pvalue_histogram.pdf")
ggsave(output_file, null_pvalue_hist, width=7.2, height=4.0, units="in")




#####################
# Figure 4: Histogram of p-values (non-null simulations)
#####################
pvalue_hist = make_pvalue_histogram(sim_df, "est_delta_beta_var_p", 0.05, "p-value", "Number of simulated genes")
output_file = paste0(simulation_visualization_dir, "no_ld_simulation_", sample_size, "_delta_beta_var_pvalue_histogram.pdf")
ggsave(output_file, pvalue_hist, width=7.2, height=4.0, units="in")




#####################
# Figure 5: Power at p-value thresholds (non-null simulations)
#####################
pvalue_thresholds = c(0.1, 0.05, 0.01)
power_summary_df = compute_power_summary_df(sim_df, "est_delta_beta_var_p", pvalue_thresholds)
power_plot = make_power_bar_plot(power_summary_df, "p-value threshold", "Power")
output_file = paste0(simulation_visualization_dir, "no_ld_simulation_", sample_size, "_delta_beta_var_power.pdf")
ggsave(output_file, power_plot, width=7.2, height=4.0, units="in")
