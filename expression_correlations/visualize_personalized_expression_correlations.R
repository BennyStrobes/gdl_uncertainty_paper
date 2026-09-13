args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(scales)
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


# Expression-FSR bins shared by the calibration and mean correlation plots
fsr_bin_breaks = c(0, 0.1, 0.2, 0.3, 0.4, Inf)
fsr_bin_labels = c("[0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", ">0.4")

# Heritable genes: cis-SNP heritability LRT p-value below this threshold
heritability_pvalue_threshold = 0.05


load_per_tissue_expression_correlations <- function(per_tissue_expression_correlation_dir, tissue_info_df) {
	# Stack every tissue's per-gene output file
	# (<tissue>_<sample>_personalized_expression_prediction.txt) into one long data frame with a
	# target_tissue column. Tissues without a results file are skipped with a warning.
	results_df = data.frame()
	for (row_iter in seq_len(nrow(tissue_info_df))) {
		target_tissue = tissue_info_df$GTEx_tissue[row_iter]
		target_sample = tissue_info_df$target_identifier[row_iter]
		results_file = paste0(per_tissue_expression_correlation_dir, target_tissue, "_", target_sample, "_personalized_expression_prediction.txt")
		if (file.exists(results_file) == FALSE) {
			print(paste("WARNING: per-tissue expression correlation file not found; skipping:", results_file))
			next
		}
		tissue_df = read.table(results_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tissue_df$target_tissue = target_tissue
		results_df = rbind(results_df, tissue_df)
	}
	if (nrow(results_df) == 0) {
		stop("No per-tissue expression correlation files found")
	}
	return(results_df)
}


restrict_to_heritable_genes <- function(df, pvalue_threshold) {
	# Keep genes whose cis-SNP heritability is significant at pvalue_threshold
	heritable_df = df[!is.na(df$cis_snp_h2_pvalue) & df$cis_snp_h2_pvalue < pvalue_threshold, ]
	print(paste("Restricting to heritable genes (p <", pvalue_threshold, "):", nrow(heritable_df), "of", nrow(df), "genes kept"))
	return(heritable_df)
}


compute_fsr_bin_summary_df <- function(df, fsr_col, tissue_colors) {
	# One row per (tissue, expression-FSR bin): observed FSR (fraction of genes whose rescaled
	# expression correlation is negative) with a binomial 95% CI, expected FSR (mean of the
	# per-gene FSR in the bin), and the mean rescaled expression correlation with a 95% CI.
	df = df[df$target_tissue %in% names(tissue_colors), ]
	df = df[!is.na(df[[fsr_col]]) & !is.na(df$rescaled_expression_correlation), ]
	df$fsr_bin = cut(df[[fsr_col]], breaks=fsr_bin_breaks, labels=fsr_bin_labels, include.lowest=TRUE, right=TRUE)

	tissue_arr = c()
	bin_arr = c()
	x_arr = c()
	n_gene_arr = c()
	observed_fsr_arr = c()
	observed_fsr_lb_arr = c()
	observed_fsr_ub_arr = c()
	expected_fsr_arr = c()
	mean_corr_arr = c()
	mean_corr_lb_arr = c()
	mean_corr_ub_arr = c()

	for (target_tissue in names(tissue_colors)) {
		for (bin_iter in 1:length(fsr_bin_labels)) {
			bin_name = fsr_bin_labels[bin_iter]
			bin_df = df[df$target_tissue == target_tissue & df$fsr_bin == bin_name, ]
			n_genes = nrow(bin_df)

			if (n_genes == 0) {
				observed_fsr = NA
				observed_fsr_lb = NA
				observed_fsr_ub = NA
				expected_fsr = NA
				mean_corr = NA
				mean_corr_lb = NA
				mean_corr_ub = NA
			} else {
				n_false_sign = sum(bin_df$rescaled_expression_correlation < 0)
				observed_fsr = n_false_sign/n_genes
				observed_fsr_se = sqrt(observed_fsr*(1.0 - observed_fsr)/n_genes)
				observed_fsr_lb = max(0.0, observed_fsr - 1.96*observed_fsr_se)
				observed_fsr_ub = min(1.0, observed_fsr + 1.96*observed_fsr_se)
				expected_fsr = mean(bin_df[[fsr_col]])
				mean_corr = mean(bin_df$rescaled_expression_correlation)
				if (n_genes > 1) {
					mean_corr_se = sd(bin_df$rescaled_expression_correlation)/sqrt(n_genes)
				} else {
					mean_corr_se = NA
				}
				mean_corr_lb = mean_corr - 1.96*mean_corr_se
				mean_corr_ub = mean_corr + 1.96*mean_corr_se
			}

			tissue_arr = c(tissue_arr, target_tissue)
			bin_arr = c(bin_arr, bin_name)
			x_arr = c(x_arr, bin_iter)
			n_gene_arr = c(n_gene_arr, n_genes)
			observed_fsr_arr = c(observed_fsr_arr, observed_fsr)
			observed_fsr_lb_arr = c(observed_fsr_lb_arr, observed_fsr_lb)
			observed_fsr_ub_arr = c(observed_fsr_ub_arr, observed_fsr_ub)
			expected_fsr_arr = c(expected_fsr_arr, expected_fsr)
			mean_corr_arr = c(mean_corr_arr, mean_corr)
			mean_corr_lb_arr = c(mean_corr_lb_arr, mean_corr_lb)
			mean_corr_ub_arr = c(mean_corr_ub_arr, mean_corr_ub)
		}
	}

	summary_df = data.frame(
		target_tissue=tissue_arr,
		fsr_bin=bin_arr,
		x=x_arr,
		n_genes=n_gene_arr,
		observed_fsr=observed_fsr_arr,
		observed_fsr_lb=observed_fsr_lb_arr,
		observed_fsr_ub=observed_fsr_ub_arr,
		expected_fsr=expected_fsr_arr,
		mean_corr=mean_corr_arr,
		mean_corr_lb=mean_corr_lb_arr,
		mean_corr_ub=mean_corr_ub_arr
	)
	summary_df$fsr_bin = factor(summary_df$fsr_bin, levels=fsr_bin_labels)
	tissue_levels = names(tissue_colors)
	summary_df$tissue = factor(summary_df$target_tissue, levels=tissue_levels, labels=gsub("_", " ", tissue_levels))
	return(summary_df)
}


make_per_tissue_fsr_calibration_panel_plot <- function(summary_df, tissue_colors, xlab) {
	# Calibration plot with one panel per tissue (side by side): observed FSR (points with binomial
	# 95% CI, joined by a line, colored by tissue) against the expected FSR (dashed segment) in each
	# expression-FSR bin. Each panel's x-axis labels carry that tissue's per-bin gene counts.
	reference_color = "#3F3F46"
	tissue_levels = levels(summary_df$tissue)
	tissue_colors_use = as.character(tissue_colors[names(tissue_colors)])

	# Per-panel x-axis labels (bin + N) through a labeller-free trick: facet on tissue, and draw the
	# gene counts as text under the points instead of in the tick labels (tick labels must be shared).
	summary_df$n_label = paste0("N=", summary_df$n_genes)

	pp = ggplot(summary_df, aes(x=x, y=observed_fsr, color=tissue)) +
		geom_segment(aes(x=x - 0.35, xend=x + 0.35, y=expected_fsr, yend=expected_fsr), color=reference_color, linetype="dashed", linewidth=0.6, na.rm=TRUE) +
		geom_line(linewidth=0.7, na.rm=TRUE) +
		geom_pointrange(aes(ymin=observed_fsr_lb, ymax=observed_fsr_ub), linewidth=0.7, na.rm=TRUE) +
		geom_point(size=2.4, na.rm=TRUE) +
		geom_text(aes(y=-0.06, label=n_label), color=reference_color, size=2.4, na.rm=TRUE) +
		facet_wrap(~tissue, nrow=1) +
		scale_color_manual(values=tissue_colors_use) +
		scale_x_continuous(breaks=seq_along(fsr_bin_labels), labels=fsr_bin_labels) +
		scale_y_continuous(limits=c(-0.09, 0.6), breaks=seq(0, 0.6, by=0.1)) +
		figure_theme() +
		theme(
			legend.position="none",
			strip.background=element_blank(),
			strip.text=element_text(face="bold", size=11),
			axis.text.x=element_text(angle=35, hjust=1, vjust=1, size=8),
			panel.spacing=unit(0.8, "lines")
		) +
		labs(x=xlab, y="Observed FSR")
	return(pp)
}


make_overlaid_fsr_calibration_plot <- function(summary_df, tissue_colors, xlab) {
	# Calibration plot with all tissues in one panel: within each expression-FSR bin the tissues are
	# dodged, each with its observed FSR (point + binomial 95% CI, colored by tissue) and its own
	# expected FSR drawn as a short dashed segment in the same color.
	tissue_levels = levels(summary_df$tissue)
	tissue_colors_use = as.character(tissue_colors[names(tissue_colors)])
	n_tissues = length(tissue_levels)
	dodge_width = 0.7
	# Manual dodge so the expected-FSR segments line up with the dodged points
	tissue_offsets = seq(-dodge_width/2, dodge_width/2, length.out=n_tissues)
	summary_df$x_dodged = summary_df$x + tissue_offsets[as.integer(summary_df$tissue)]
	segment_half_width = (dodge_width/(n_tissues - 1))*0.4

	pp = ggplot(summary_df, aes(x=x_dodged, y=observed_fsr, color=tissue)) +
		geom_segment(aes(x=x_dodged - segment_half_width, xend=x_dodged + segment_half_width, y=expected_fsr, yend=expected_fsr), linetype="dashed", linewidth=0.6, na.rm=TRUE) +
		geom_line(aes(group=tissue), linewidth=0.5, alpha=0.6, na.rm=TRUE) +
		geom_pointrange(aes(ymin=observed_fsr_lb, ymax=observed_fsr_ub), linewidth=0.6, na.rm=TRUE) +
		geom_point(size=2.2, na.rm=TRUE) +
		scale_color_manual(values=tissue_colors_use, name="Tissue") +
		scale_x_continuous(breaks=seq_along(fsr_bin_labels), labels=fsr_bin_labels) +
		scale_y_continuous(limits=c(0, 0.6), breaks=seq(0, 0.6, by=0.1)) +
		figure_theme() +
		theme(
			legend.position="right",
			legend.title=element_text(face="bold"),
			plot.margin=margin(8, 14, 8, 8)
		) +
		labs(x=xlab, y="Observed FSR")
	return(pp)
}


make_five_tissue_mean_correlation_bar_plot <- function(summary_df, tissue_colors, xlab) {
	# Grouped bar plot: x-axis is the expression-FSR bin and within each bin one dodged bar per tissue
	# giving the mean rescaled expression correlation across genes (95% CI error bars).
	tissue_colors_use = as.character(tissue_colors[names(tissue_colors)])
	pp = ggplot(summary_df, aes(x=fsr_bin, y=mean_corr, fill=tissue)) +
		geom_col(position=position_dodge(width=.8), width=.72, color="#111827", linewidth=.2, na.rm=TRUE) +
		geom_errorbar(aes(ymin=mean_corr_lb, ymax=mean_corr_ub), position=position_dodge(width=.8), width=.18, linewidth=.35, color="#111827", na.rm=TRUE) +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_fill_manual(values=tissue_colors_use, name="Tissue") +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		figure_theme() +
		theme(
			legend.position="right",
			legend.title=element_text(face="bold"),
			plot.margin=margin(8, 14, 8, 8)
		) +
		labs(x=xlab, y="Average expression\ncorrelation")
	return(pp)
}


make_stacked_shared_x_plot <- function(top_plot, bottom_plot, panel_labels) {
	# Stack two plots that share the same x-axis into one figure with a single legend at the bottom:
	# the top panel drops its (duplicated) x-axis title, tick labels and ticks.
	panel_tag_theme = theme(plot.tag=element_text(face="bold", size=14, hjust=0, vjust=1), plot.tag.position=c(0, 1))
	top_panel = top_plot + labs(tag=panel_labels[1]) + panel_tag_theme + theme(
		legend.position="none",
		axis.title.x=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		plot.margin=margin(6, 8, 2, 8)
	)
	bottom_panel = bottom_plot + labs(tag=panel_labels[2]) + panel_tag_theme + theme(legend.position="none", plot.margin=margin(2, 8, 4, 8))
	top_grob = ggplotGrob(top_panel)
	bottom_grob = ggplotGrob(bottom_panel)
	shared_widths = grid::unit.pmax(top_grob$widths, bottom_grob$widths)
	top_grob$widths = shared_widths
	bottom_grob$widths = shared_widths
	stacked_panels = ggdraw(rbind(top_grob, bottom_grob, size="first"))

	legend_source_plot = bottom_plot +
		theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=9), legend.margin=margin(0, 0, 0, 0)) +
		guides(fill=guide_legend(nrow=1))
	# ggplot2 >= 3.5 names the bottom guide box "guide-box-bottom"; older versions have a single "guide-box"
	shared_legend_grob = tryCatch(get_plot_component(legend_source_plot, "guide-box-bottom"), error=function(e) NULL)
	if (is.null(shared_legend_grob) || inherits(shared_legend_grob, "zeroGrob")) {
		shared_legend_grob = get_plot_component(legend_source_plot, "guide-box")
	}
	return(plot_grid(stacked_panels, shared_legend_grob, ncol=1, rel_heights=c(1, .08)))
}





#####################
# Command line args
#####################
per_tissue_expression_correlation_dir = args[1]
tissue_names_file = args[2]
visualization_dir = args[3]


#####################
# Load in data
#####################
# Load in tissue info df
tissue_info_df = read.table(tissue_names_file, header=TRUE, sep="\t")

# Per-gene expression correlation results across all tissues
results_df = load_per_tissue_expression_correlations(per_tissue_expression_correlation_dir, tissue_info_df)

# Restrict to heritable genes
heritable_results_df = restrict_to_heritable_genes(results_df, heritability_pvalue_threshold)


#####################
# Five tissues used for plotting (same colors as the S-LDMC five-tissue plots)
# Heart left ventricle: #4C72B0
# Brain Cortex #DD8452
# Liver #55A868
# Whole Blood #C44E52
# Muscle Skeletal #8172B3
#####################
five_tissue_colors = c(
	"Heart_Left_Ventricle"="#4C72B0",
	"Brain_Cortex"="#DD8452",
	"Liver"="#55A868",
	"Whole_Blood"="#C44E52",
	"Muscle_Skeletal"="#8172B3"
)
missing_tissues = setdiff(names(five_tissue_colors), unique(heritable_results_df$target_tissue))
if (length(missing_tissues) > 0) {
	print(paste("WARNING: no expression correlation results for:", paste(missing_tissues, collapse=", ")))
}


#####################
# Calibration and mean correlation plots, once per expression-FSR definition
# expression_FSR: residual variance constant within borzoi magnitude bin
# expression_FSR_af_specific: allele-frequency-specific residual variance
#####################
fsr_definitions = c("expression_FSR", "expression_FSR_af_specific")
fsr_xlabs = c("expression_FSR"="Expression-FSR", "expression_FSR_af_specific"="Expression-FSR (AF-specific)")

for (fsr_col in fsr_definitions) {
	print(paste("Plotting", fsr_col))
	summary_df = compute_fsr_bin_summary_df(heritable_results_df, fsr_col, five_tissue_colors)
	print(summary_df)
	xlab = fsr_xlabs[[fsr_col]]

	# Calibration: one panel per tissue
	per_tissue_calibration_plot = make_per_tissue_fsr_calibration_panel_plot(summary_df, five_tissue_colors, xlab)
	ggsave(paste0(visualization_dir, "five_tissue_", fsr_col, "_calibration_per_tissue_panels.pdf"), per_tissue_calibration_plot, width=9.5, height=2.9)

	# Calibration: all tissues overlaid in one panel
	overlaid_calibration_plot = make_overlaid_fsr_calibration_plot(summary_df, five_tissue_colors, xlab)
	ggsave(paste0(visualization_dir, "five_tissue_", fsr_col, "_calibration_overlaid.pdf"), overlaid_calibration_plot, width=6.0, height=3.2)

	# Mean rescaled expression correlation in each expression-FSR bin
	mean_correlation_plot = make_five_tissue_mean_correlation_bar_plot(summary_df, five_tissue_colors, xlab)
	ggsave(paste0(visualization_dir, "five_tissue_", fsr_col, "_mean_correlation.pdf"), mean_correlation_plot, width=6.0, height=3.2)

	# Joint calibration (overlaid) + mean correlation plot (shared x-axis and shared legend)
	joint_plot = make_stacked_shared_x_plot(overlaid_calibration_plot, mean_correlation_plot, c("a", "b"))
	ggsave(paste0(visualization_dir, "five_tissue_", fsr_col, "_joint_calibration_mean_correlation.pdf"), joint_plot, width=6.0, height=5.4)
}
