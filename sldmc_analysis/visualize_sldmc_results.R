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



make_sldmc_magnitude_bar_plot <- function(df, output_name, ylab, plot_title, bar_color) {
	# Bar plot of one output variable (correlation, calibration_slope, ...) across the borzoi
	# magnitude bins. Expects the newer bootstrap_stats format with annotation_name/category_name,
	# so the magnitude bins live in the category_name column of the borzoi_magnitude_bins annotation.
	plot_df = df[df$annotation_name == "borzoi_magnitude_bins" & df$output_name == output_name, ]
	if (nrow(plot_df) == 0) {
		stop(paste("No borzoi_magnitude_bins rows found for output", output_name))
	}
	bin_ids = suppressWarnings(as.numeric(sub(".*magnitude_bin([0-9]+)$", "\\1", plot_df$category_name)))
	plot_df = plot_df[order(bin_ids), ]
	plot_df$category_name = factor(plot_df$category_name, levels=plot_df$category_name)
	magnitude_bin_labels = c("<0.001", "0.001-0.01", "0.01-0.075", "0.075-0.2", ">=0.2")
	if (nrow(plot_df) == length(magnitude_bin_labels)) {
		plot_df$annotation_label = magnitude_bin_labels
	} else {
		plot_df$annotation_label = paste0("bin", seq_len(nrow(plot_df)))
	}
	# Gaussian approximation confidence intervals
	plot_df$ci_lower = plot_df$mean - 1.96*plot_df$bootstrap_se
	plot_df$ci_upper = plot_df$mean + 1.96*plot_df$bootstrap_se
	return(
		ggplot(plot_df, aes(x=category_name, y=mean)) +
		geom_col(fill=bar_color, color="#111827", width=.72, linewidth=.25) +
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.18, linewidth=.45, color="#111827") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		scale_x_discrete(labels=plot_df$annotation_label) +
		labs(title=plot_title) +
		xlab("Borzoi magnitude bin") +
		ylab(ylab) +
		figure_theme() +
		theme(
			axis.text.x=element_text(angle=35, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}


make_sldmc_bin_bar_plot <- function(df, annotation_name, output_name, ylab, plot_title, xlab_title, bar_color, bin_labels, y_accuracy=.01) {
	# Bar plot of one output variable across the bins of a single binned annotation (e.g. af_bins,
	# dist_to_tss_bins). The bins live in the category_name column; bin ordering comes from the
	# trailing bin index, and bin_labels (if the count matches) provide human-readable x-axis labels.
	plot_df = df[df$annotation_name == annotation_name & df$output_name == output_name, ]
	if (nrow(plot_df) == 0) {
		stop(paste("No", annotation_name, "rows found for output", output_name))
	}
	bin_ids = suppressWarnings(as.numeric(sub(".*[^0-9]([0-9]+)$", "\\1", plot_df$category_name)))
	plot_df = plot_df[order(bin_ids), ]
	plot_df$category_name = factor(plot_df$category_name, levels=plot_df$category_name)
	if (nrow(plot_df) == length(bin_labels)) {
		plot_df$annotation_label = bin_labels
	} else {
		plot_df$annotation_label = paste0("bin", seq_len(nrow(plot_df)))
	}
	# Gaussian approximation confidence intervals
	plot_df$ci_lower = plot_df$mean - 1.96*plot_df$bootstrap_se
	plot_df$ci_upper = plot_df$mean + 1.96*plot_df$bootstrap_se
	return(
		ggplot(plot_df, aes(x=category_name, y=mean)) +
		geom_col(fill=bar_color, color="#111827", width=.72, linewidth=.25) +
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.18, linewidth=.45, color="#111827") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_y_continuous(labels=number_format(accuracy=y_accuracy)) +
		scale_x_discrete(labels=plot_df$annotation_label) +
		labs(title=plot_title) +
		xlab(xlab_title) +
		ylab(ylab) +
		figure_theme() +
		theme(
			axis.text.x=element_text(angle=35, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}


make_sldmc_magnitude_stratified_bar_plot <- function(df, base_annotation_name, output_name, ylab, base_legend_title, bar_color, base_bin_labels) {
	# Grouped bar plot: x-axis is the borzoi magnitude bin, and within each magnitude bin the base
	# annotation's bins (e.g. distance to TSS) are shown as separate colored, dodged bars (darker =
	# higher base bin). Reads the magnitude_stratified file: annotation borzoi_magnitude_stratifiedX
	# <base>, categories magnitude_bin<m>X<base_bin>.
	stratified_annotation_name = paste0("borzoi_magnitude_stratifiedX", base_annotation_name)
	plot_df = df[df$annotation_name == stratified_annotation_name & df$output_name == output_name, ]
	if (nrow(plot_df) == 0) {
		stop(paste("No", stratified_annotation_name, "rows found for output", output_name))
	}
	plot_df$magnitude_bin_id = suppressWarnings(as.numeric(sub("^magnitude_bin([0-9]+)X.*", "\\1", plot_df$category_name)))
	plot_df$base_bin_id = suppressWarnings(as.numeric(sub(".*[^0-9]([0-9]+)$", "\\1", plot_df$category_name)))

	magnitude_bin_labels = c("<0.001", "0.001-0.01", "0.01-0.075", "0.075-0.2", ">=0.2")
	magnitude_ids = sort(unique(plot_df$magnitude_bin_id))
	base_ids = sort(unique(plot_df$base_bin_id))
	magnitude_labels_use = if (length(magnitude_ids) == length(magnitude_bin_labels)) magnitude_bin_labels else paste0("mag bin", magnitude_ids)
	base_labels_use = if (length(base_ids) == length(base_bin_labels)) base_bin_labels else paste0("bin", base_ids)

	plot_df$magnitude_bin = factor(plot_df$magnitude_bin_id, levels=magnitude_ids, labels=magnitude_labels_use)
	plot_df$base_bin = factor(plot_df$base_bin_id, levels=base_ids, labels=base_labels_use)
	# Gaussian approximation confidence intervals
	plot_df$ci_lower = plot_df$mean - 1.96*plot_df$bootstrap_se
	plot_df$ci_upper = plot_df$mean + 1.96*plot_df$bootstrap_se

	# Shades of the metric's base color, one per base-annotation bin (light -> bar_color)
	base_colors = colorRampPalette(c("white", bar_color))(length(base_ids) + 1)[-1]
	return(
		ggplot(plot_df, aes(x=magnitude_bin, y=mean, fill=base_bin)) +
		geom_col(position=position_dodge(width=.8), width=.72, color="#111827", linewidth=.2) +
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position=position_dodge(width=.8), width=.18, linewidth=.35, color="#111827") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_fill_manual(values=base_colors, name=base_legend_title) +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		xlab("Borzoi magnitude bin") +
		ylab(ylab) +
		figure_theme() +
		theme(
			legend.position="right",
			legend.title=element_text(face="bold"),
			axis.text.x=element_text(angle=35, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}


per_tissue_sldmc_magnitude_results <- function(sldmc_results_output_dir, tissue_info_df, annotation_version="default") {
	# Build a long data frame with one row per (tissue, borzoi magnitude bin), holding both the
	# correlation and the calibration slope (with bootstrap SEs) from that tissue's per-tissue
	# S-LDMC bootstrap_stats file (sldmc_results_<tissue>_<sample>_<version>_bootstrap_stats.txt).
	results_df = data.frame()
	for (row_iter in seq_len(nrow(tissue_info_df))) {
		target_tissue = tissue_info_df$GTEx_tissue[row_iter]
		target_sample = tissue_info_df$target_identifier[row_iter]
		stats_file = paste0(
			sldmc_results_output_dir,
			"sldmc_results_",
			target_tissue,
			"_",
			target_sample,
			"_",
			annotation_version,
			"_bootstrap_stats.txt"
		)
		if (file.exists(stats_file) == FALSE) {
			print(paste("WARNING: per-tissue S-LDMC stats file not found; skipping:", stats_file))
			next
		}
		tissue_df = read.table(stats_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		magnitude_df = tissue_df[tissue_df$annotation_name == "borzoi_magnitude_bins", ]
		if (nrow(magnitude_df) == 0) {
			print(paste("WARNING: no borzoi_magnitude_bins rows in", stats_file))
			next
		}

		correlation_df = magnitude_df[magnitude_df$output_name == "correlation", c("category_name", "mean", "bootstrap_se")]
		colnames(correlation_df) = c("category_name", "correlation", "correlation_se")
		calibration_df = magnitude_df[magnitude_df$output_name == "calibration_slope", c("category_name", "mean", "bootstrap_se")]
		colnames(calibration_df) = c("category_name", "calibration_slope", "calibration_slope_se")

		merged_df = merge(correlation_df, calibration_df, by="category_name")
		bin_ids = suppressWarnings(as.numeric(sub(".*magnitude_bin([0-9]+)$", "\\1", merged_df$category_name)))
		merged_df = merged_df[order(bin_ids), ]
		merged_df$target_tissue = target_tissue
		merged_df$target_sample = target_sample
		merged_df = merged_df[, c("target_tissue", "target_sample", "category_name", "correlation", "correlation_se", "calibration_slope", "calibration_slope_se")]
		results_df = rbind(results_df, merged_df)
	}
	if (nrow(results_df) == 0) {
		stop(paste("No per-tissue S-LDMC magnitude results loaded from", sldmc_results_output_dir))
	}
	return(results_df)
}


make_per_tissue_magnitude_dot_plot <- function(df, value_col, ylab, point_color) {
	# Jittered dot plot of a per-tissue S-LDMC estimate (correlation, calibration_slope, ...) across
	# the borzoi magnitude bins (one dot per tissue per bin). Expects the long df from
	# per_tissue_sldmc_magnitude_results.
	plot_df = df
	plot_df$value = plot_df[[value_col]]
	bin_ids = suppressWarnings(as.numeric(sub(".*magnitude_bin([0-9]+)$", "\\1", plot_df$category_name)))
	ordered_bins = unique(plot_df$category_name[order(bin_ids)])
	plot_df$category_name = factor(plot_df$category_name, levels=ordered_bins)
	magnitude_bin_labels = c("<0.001", "0.001-0.01", "0.01-0.075", "0.075-0.2", ">=0.2")
	if (nlevels(plot_df$category_name) == length(magnitude_bin_labels)) {
		bin_labels = magnitude_bin_labels
	} else {
		bin_labels = paste0("bin", seq_len(nlevels(plot_df$category_name)))
	}
	return(
		ggplot(plot_df, aes(x=category_name, y=value)) +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		geom_jitter(position=position_jitter(width=.15, height=0, seed=1), shape=21, size=1.8, stroke=.35, fill=point_color, color="#111827", alpha=.8) +
		scale_y_continuous(labels=number_format(accuracy=.01)) +
		scale_x_discrete(labels=bin_labels) +
		xlab("Borzoi magnitude bin") +
		ylab(ylab) +
		figure_theme() +
		theme(
			axis.text.x=element_text(angle=35, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}


make_five_tissue_magnitude_bar_plot <- function(df, value_col, se_col, ylab, tissue_colors, y_limits=NULL, y_accuracy=.01, reference_line=NULL, truncation_note=NULL) {
	# Grouped bar plot: x-axis is the borzoi magnitude bin, and within each magnitude bin one dodged
	# bar per selected tissue (colored by tissue, with 95% Gaussian CI error bars). Expects the long df
	# from per_tissue_sldmc_magnitude_results and a named vector mapping target_tissue -> bar color;
	# the order of that vector sets the bar / legend order.
	# y_limits zooms the y-axis without dropping data (error bars running past it are drawn to the panel
	# edge), y_accuracy sets the tick label precision, reference_line adds a dotted horizontal reference
	# (e.g. 1 = perfect calibration), and truncation_note labels the panel when y_limits clips a CI.
	plot_df = df[df$target_tissue %in% names(tissue_colors), ]
	missing_tissues = setdiff(names(tissue_colors), unique(plot_df$target_tissue))
	if (length(missing_tissues) > 0) {
		print(paste("WARNING: no per-tissue S-LDMC results for:", paste(missing_tissues, collapse=", ")))
	}
	if (nrow(plot_df) == 0) {
		stop("No per-tissue S-LDMC rows found for any of the selected tissues")
	}
	plot_df$value = plot_df[[value_col]]
	plot_df$value_se = plot_df[[se_col]]

	bin_ids = suppressWarnings(as.numeric(sub(".*magnitude_bin([0-9]+)$", "\\1", plot_df$category_name)))
	ordered_bins = unique(plot_df$category_name[order(bin_ids)])
	plot_df$category_name = factor(plot_df$category_name, levels=ordered_bins)
	magnitude_bin_labels = c("<0.001", "0.001-0.01", "0.01-0.075", "0.075-0.2", ">=0.2")
	if (nlevels(plot_df$category_name) == length(magnitude_bin_labels)) {
		bin_labels = magnitude_bin_labels
	} else {
		bin_labels = paste0("bin", seq_len(nlevels(plot_df$category_name)))
	}

	tissue_levels = names(tissue_colors)[names(tissue_colors) %in% plot_df$target_tissue]
	plot_df$tissue = factor(plot_df$target_tissue, levels=tissue_levels, labels=gsub("_", " ", tissue_levels))
	tissue_colors_use = as.character(tissue_colors[tissue_levels])
	# Gaussian approximation confidence intervals
	plot_df$ci_lower = plot_df$value - 1.96*plot_df$value_se
	plot_df$ci_upper = plot_df$value + 1.96*plot_df$value_se

	# Optional layers (a NULL layer is a no-op when added to a ggplot)
	reference_line_layer = NULL
	if (is.null(reference_line) == FALSE) {
		reference_line_layer = geom_hline(yintercept=reference_line, linewidth=.4, color="#6B7280", linetype="dotted")
	}
	truncation_note_layer = NULL
	if (is.null(truncation_note) == FALSE) {
		truncation_note_layer = annotate("text", x=Inf, y=Inf, label=truncation_note, hjust=1.04, vjust=1.6, size=2.9, color="#6B7280")
	}
	coord_layer = NULL
	if (is.null(y_limits) == FALSE) {
		coord_layer = coord_cartesian(ylim=y_limits)
	}

	return(
		ggplot(plot_df, aes(x=category_name, y=value, fill=tissue)) +
		geom_col(position=position_dodge(width=.8), width=.72, color="#111827", linewidth=.2) +
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position=position_dodge(width=.8), width=.18, linewidth=.35, color="#111827") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		reference_line_layer +
		truncation_note_layer +
		coord_layer +
		scale_fill_manual(values=tissue_colors_use, name="Tissue") +
		scale_y_continuous(labels=number_format(accuracy=y_accuracy)) +
		scale_x_discrete(labels=bin_labels) +
		xlab("Borzoi magnitude bin") +
		ylab(ylab) +
		figure_theme() +
		theme(
			legend.position="right",
			legend.title=element_text(face="bold"),
			axis.text.x=element_text(angle=35, hjust=1, vjust=1),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}


make_stacked_shared_x_plot <- function(top_plot, bottom_plot, panel_labels, shared_legend=TRUE, horizontal_x_labels=TRUE) {
	# Stack two plots that share the same x-axis into one figure: the top panel drops its (duplicated)
	# x-axis title, tick labels and ticks, so the axis is drawn once under the bottom panel.
	# shared_legend replaces the two per-panel legends with a single legend along the bottom (set FALSE
	# for plots that have no legend), and horizontal_x_labels un-angles the bottom panel's tick labels
	# (set FALSE when the labels are too long to fit horizontally).
	# The panel letters are ggplot tags rather than plot_grid labels so that they survive the gtable
	# stacking below.
	panel_tag_theme = theme(plot.tag=element_text(face="bold", size=14, hjust=0, vjust=1), plot.tag.position=c(0, 1))
	top_panel = top_plot + labs(tag=panel_labels[1]) + panel_tag_theme + theme(
		legend.position="none",
		axis.title.x=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		plot.margin=margin(6, 8, 2, 8)
	)
	bottom_panel = bottom_plot + labs(tag=panel_labels[2]) + panel_tag_theme + theme(legend.position="none", plot.margin=margin(2, 8, 4, 8))
	if (horizontal_x_labels) {
		bottom_panel = bottom_panel + theme(axis.text.x=element_text(angle=0, hjust=.5, vjust=1))
	}
	# Stack the two plots as gtables instead of via plot_grid: plot_grid gives each plot an equal-height
	# cell, so the panel that still carries the x-axis labels ends up with a shorter plotting area. In a
	# combined gtable the two panel rows are the only "null" heights, so they split the leftover space
	# evenly and both panels get an identical y-axis height.
	top_grob = ggplotGrob(top_panel)
	bottom_grob = ggplotGrob(bottom_panel)
	shared_widths = grid::unit.pmax(top_grob$widths, bottom_grob$widths)
	top_grob$widths = shared_widths
	bottom_grob$widths = shared_widths
	stacked_panels = ggdraw(rbind(top_grob, bottom_grob, size="first"))
	if (shared_legend == FALSE) {
		return(stacked_panels)
	}
	legend_source_plot = bottom_plot +
		theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=9), legend.margin=margin(0, 0, 0, 0)) +
		guides(fill=guide_legend(nrow=1))
	# ggplot2 >= 3.5 names the bottom guide box "guide-box-bottom"; older versions have a single "guide-box"
	shared_legend_grob = tryCatch(get_plot_component(legend_source_plot, "guide-box-bottom"), error=function(e) NULL)
	if (is.null(shared_legend_grob) || inherits(shared_legend_grob, "zeroGrob")) {
		shared_legend_grob = get_plot_component(legend_source_plot, "guide-box")
	}
	return(plot_grid(stacked_panels, shared_legend_grob, nrow=2, rel_heights=c(1, .09)))
}


subset_diff_by_source <- function(diff_df, anno_to_source, target_source) {
	# Keep only the diff-stats rows whose annotation is of a given source (sldsc / gene_set). The
	# diff file has no source column, so source is looked up by (base) annotation name in anno_to_source.
	# Magnitude-stratified names carry a prefix, so strip it before the lookup.
	base_names = sub("^borzoi_magnitude_stratifiedX", "", diff_df$annotation_name)
	row_sources = anno_to_source[base_names]
	keep = !is.na(row_sources) & row_sources == target_source
	return(diff_df[keep, ])
}


save_source_forest_plot <- function(diff_df, anno_to_source, target_source, excluded_annotation_regex, present_category, output_file) {
	# Build and save a forest plot for a single annotation source, with its own (independent) FDR.
	subset_df = subset_diff_by_source(diff_df, anno_to_source, target_source)
	present_rows = subset_df$category_name == present_category
	base_names = sub("^borzoi_magnitude_stratifiedX", "", subset_df$annotation_name)
	excluded = grepl(excluded_annotation_regex, base_names, ignore.case=TRUE)
	n_annotations = length(unique(subset_df$annotation_name[present_rows & (excluded == FALSE)]))
	if (n_annotations == 0) {
		print(paste("WARNING: no", target_source, "annotations to plot; skipping", output_file))
		return(invisible(NULL))
	}
	forest_plot = make_sldsc_annotation_forest_plot(subset_df, excluded_annotation_regex, present_category)
	ggsave(output_file, forest_plot, width=7.6, height=max(4, 0.22*n_annotations + 1.4))
	print(output_file)
}


make_sldsc_annotation_forest_plot <- function(diff_df, excluded_annotation_regex, present_category="present") {
	# Two-panel forest plot of the S-LDSC (binary) annotations, showing each annotation's difference
	# from the intercept for correlation (left) and calibration slope (right). Binary annotations are
	# the ones with a "present" cell; flanking / allele-frequency annotations are dropped via
	# excluded_annotation_regex. present_category selects which "present" cell to plot: "present" for
	# the default file, or "magnitude_bin<m>Xpresent" for the magnitude-stratified file (i.e. the
	# difference conditional on being in magnitude bin <m>). Points are colored by BH-FDR significance.
	present_df = diff_df[diff_df$category_name == present_category & diff_df$output_name %in% c("correlation", "calibration_slope"), ]
	# Exclusions match on the base annotation name (magnitude-stratified names carry a prefix), so
	# anchored patterns like "^MohammadiVg" work in both the default and magnitude-stratified files.
	present_base_names = sub("^borzoi_magnitude_stratifiedX", "", present_df$annotation_name)
	present_df = present_df[grepl(excluded_annotation_regex, present_base_names, ignore.case=TRUE) == FALSE, ]
	if (nrow(present_df) == 0) {
		stop("No annotation present-cell rows found for the forest plot")
	}

	# Gaussian CI on the difference, and BH-FDR significance within each metric
	present_df$ci_lower = present_df$difference - 1.96*present_df$difference_se
	present_df$ci_upper = present_df$difference + 1.96*present_df$difference_se
	present_df$pvalue = 2*pnorm(-abs(present_df$gaussian_z_score))
	present_df$significant = FALSE
	for (metric_name in unique(present_df$output_name)) {
		metric_rows = present_df$output_name == metric_name
		present_df$significant[metric_rows] = p.adjust(present_df$pvalue[metric_rows], method="BH") < 0.05
	}

	# Order annotations by their correlation difference (shared across both panels)
	correlation_df = present_df[present_df$output_name == "correlation", ]
	ordered_annotations = unique(c(correlation_df$annotation_name[order(correlation_df$difference)], present_df$annotation_name))
	present_df$annotation_name = factor(present_df$annotation_name, levels=ordered_annotations)
	present_df$metric = factor(
		present_df$output_name,
		levels=c("correlation", "calibration_slope"),
		labels=c("Correlation", "Calibration coefficient")
	)

	return(
		ggplot(present_df, aes(x=difference, y=annotation_name, color=significant)) +
		geom_vline(xintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper), height=0, linewidth=.45) +
		geom_point(size=1.7) +
		facet_wrap(~metric, scales="free_x") +
		scale_color_manual(values=c("FALSE"="#9CA3AF", "TRUE"="#3F88C5"), labels=c("FALSE"="n.s.", "TRUE"="FDR<0.05"), name=NULL) +
		scale_x_continuous(labels=number_format(accuracy=.01)) +
		scale_y_discrete(labels=function(annotation_names) gsub("_", " ", gsub("borzoi_magnitude_stratifiedX", "", annotation_names))) +
		xlab("Difference from intercept") +
		ylab(NULL) +
		figure_theme() +
		theme(
			legend.position="top",
			axis.text.y=element_text(size=7),
			panel.grid.major.y=element_line(color="#E5E7EB", linewidth=.3),
			strip.background=element_blank(),
			strip.text=element_text(face="bold", size=10),
			plot.margin=margin(8, 14, 8, 8)
		)
	)
}

load_simulation_summary <- function(simulation_results_dir, metric="correlation", n_sims=50, eqtl_ss="489", n_anno="6")	{
	anno_name_arr = c()
	method_arr = c()
	est_arr = c()

	for (sim_iter in 1:n_sims) {
		sim_file <- paste0(simulation_results_dir, "sim", sim_iter, "_sim_eqtl_ss_", eqtl_ss, "_", n_anno, "_anno_ld_corr_results_bootstrap_stats.txt")
		tmp_df <- read.table(sim_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		# Filter to output metric
		tmp_df <- tmp_df[(as.character(tmp_df$output_name) == metric) & (as.character(tmp_df$annotation_name) == "sim_signal_prop"), ]

		# Assumes n_anno=6: the simulation's signal-proportion grid is np.linspace(0, 1, n_anno)
		expected_category_names = c("prop_0.0", "prop_0.2", "prop_0.4", "prop_0.6", "prop_0.8", "prop_1.0")
		if (identical(tmp_df$category_name, expected_category_names) == FALSE) {
			stop(paste("Unexpected category_name values in simulation file:", sim_file))
		}

		est_arr = c(est_arr, tmp_df$mean)
		method_arr = c(method_arr, rep("S-LDMC", nrow(tmp_df)))
		anno_name_arr = c(anno_name_arr, 1:nrow(tmp_df))
	}

	for (sim_iter in 1:n_sims) {
		sim_file <- paste0(simulation_results_dir, "sim", sim_iter, "_sim_eqtl_ss_", eqtl_ss, "_", n_anno, "_anno_fm_corr_results_bootstrap_stats.txt")
		tmp_df <- read.table(sim_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		# Filter to output metric
		tmp_df <- tmp_df[(as.character(tmp_df$output_name) == metric), ]

		# Assumes n_anno=6: the simulation's signal-proportion grid is np.linspace(0, 1, n_anno)
		expected_category_names = c("anno0", "anno1", "anno2", "anno3", "anno4", "anno5")
		if (identical(tmp_df$annotation_name, expected_category_names) == FALSE) {
			stop(paste("Unexpected category_name values in simulation file:", sim_file))
		}

		est_arr = c(est_arr, tmp_df$mean)
		method_arr = c(method_arr, rep("Fine-map", nrow(tmp_df)))
		anno_name_arr = c(anno_name_arr, 1:nrow(tmp_df))
	}


	df_lined <- data.frame(
		annotation_name=anno_name_arr,
		method=method_arr,
		mean=est_arr
	)

	# Collapse to one row per (annotation_name, method) by averaging over simulations. The interval is a
	# Gaussian 95% interval around that average built from the standard error of the mean across
	# simulations (sd across simulations / sqrt(n simulations)), so it describes uncertainty in the
	# simulation average rather than the spread of individual simulations.
	group_columns = c("annotation_name", "method")
	sim_means = aggregate(df_lined["mean"], by=df_lined[group_columns], FUN=mean)
	sim_sems = aggregate(df_lined["mean"], by=df_lined[group_columns], FUN=function(estimates) sd(estimates)/sqrt(length(estimates)))
	colnames(sim_sems)[colnames(sim_sems) == "mean"] = "sem"
	df <- merge(sim_means, sim_sems, by=group_columns)
	df$"95_lb" = df$mean - 1.96*df$sem
	df$"95_ub" = df$mean + 1.96*df$sem

	return(df)
}

load_in_oracle_result <- function(simulation_oracle_results_dir, column_number, n_sims=50, n_anno="6") {
	# Load in one column of the oracle results file (1-based column_number, e.g. 3 = correlation)
	est_arr = c()
	anno_arr <- c()
	for (sim_iter in 1:n_sims) {
		sim_file <- paste0(simulation_oracle_results_dir, "sim", sim_iter, "_sim_variant_gene_annotations_", n_anno, "_true_sim_effect_summary.txt")
		tmp_df <- read.table(sim_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		est_arr = c(est_arr, tmp_df[[column_number]])
		anno_arr = c(anno_arr, 1:nrow(tmp_df))
	}

	df_lined <- data.frame(
		annotation_name=anno_arr,
		mean=est_arr
	)

	# Collapse to one row per annotation_name by averaging over simulations, with a Gaussian 95%
	# interval built from the standard error of the mean across simulations (same as in
	# load_simulation_summary, but there is only one "method" here so annotation_name is the only group)
	sim_means = aggregate(df_lined["mean"], by=df_lined["annotation_name"], FUN=mean)
	sim_sems = aggregate(df_lined["mean"], by=df_lined["annotation_name"], FUN=function(estimates) sd(estimates)/sqrt(length(estimates)))
	colnames(sim_sems)[colnames(sim_sems) == "mean"] = "sem"
	df <- merge(sim_means, sim_sems, by="annotation_name")
	df$"95_lb" = df$mean - 1.96*df$sem
	df$"95_ub" = df$mean + 1.96*df$sem

	return(df)
}


make_simulation_estimate_bar_plot <- function(summary_df, oracle_df, ylab, bar_color, y_accuracy=.01) {
	# Grouped bar plot of the simulation estimates: x-axis is the simulated annotation bin, with one
	# dodged bar per method (Fine-map / S-LDMC) and 95% CI error bars taken from the across-simulation
	# SEM columns built by load_simulation_summary. The oracle ("true simulated") value for each bin is
	# overlaid as a horizontal crossbar plus an open point, so bar-vs-crossbar shows the bias.
	# bar_color is the metric's color used by the real-data plots (blue for correlation, red for
	# calibration); S-LDMC takes it directly and Fine-map takes a grey-blended version of it, so a
	# simulation panel reads in the same hue as the real-data panel for that metric.
	bar_df = summary_df[summary_df$method %in% c("Fine-map", "S-LDMC"), ]
	if (nrow(bar_df) == 0) {
		stop("No Fine-map / S-LDMC rows found in the simulation summary data frame")
	}
	truth_df = oracle_df

	# annotation_name is the 1-based annotation bin index, used directly as the x-axis label. Bin i holds
	# the variant-gene pairs whose simulated borzoi effect is np.linspace(0, 1, n_anno)[i] true signal and
	# the rest noise (simulate_est_borzoi_effects_for_correlation_experiment.py), so the bins run from
	# pure noise up to pure signal.
	annotation_ids = sort(unique(c(bar_df$annotation_name, truth_df$annotation_name)))
	bar_df$annotation_bin = factor(bar_df$annotation_name, levels=annotation_ids)
	truth_df$annotation_bin = factor(truth_df$annotation_name, levels=annotation_ids)

	# The CI bounds live in columns whose names start with a digit, so they need quoting
	bar_df$ci_lower = bar_df$"95_lb"
	bar_df$ci_upper = bar_df$"95_ub"
	bar_df$method = factor(bar_df$method, levels=c("Fine-map", "S-LDMC"))

	# Mostly-grey version of the metric color for the Fine-map bars: near-neutral grey carrying just a
	# hint of the metric hue, so it sits beside the saturated S-LDMC bar without competing with it
	fine_map_color = colorRampPalette(c("#9CA3AF", bar_color))(10)[2]
	method_colors = c(fine_map_color, bar_color)
	names(method_colors) = c("Fine-map", "S-LDMC")

	return(
		ggplot(bar_df, aes(x=annotation_bin, y=mean, fill=method)) +
		geom_col(position=position_dodge(width=.78), width=.68, color="#111827", linewidth=.25) +
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position=position_dodge(width=.78), width=.16, linewidth=.4, color="#111827") +
		geom_errorbar(data=truth_df, aes(x=annotation_bin, ymin=mean, ymax=mean, color="True simulated"), inherit.aes=FALSE, width=.58, linewidth=.9) +
		geom_point(data=truth_df, aes(x=annotation_bin, y=mean, color="True simulated"), inherit.aes=FALSE, shape=21, size=1.8, stroke=.55, fill="white") +
		geom_hline(yintercept=0, linewidth=.4, color="#6B7280", linetype="dashed") +
		scale_fill_manual(values=method_colors) +
		scale_color_manual(values=c("True simulated"="#111827"), guide="none") +
		scale_y_continuous(labels=number_format(accuracy=y_accuracy)) +
		labs(fill=NULL, color=NULL) +
		xlab("Simulated annotation category") +
		ylab(ylab) +
		figure_theme() +
		theme(
			legend.position="top",
			legend.text=element_text(size=10),
			legend.key.size=unit(.9, "lines"),
			# Pull the legend in tight: legend.box.spacing is the gap between the panel and the legend,
			# and legend.margin the padding inside the legend box
			legend.box.spacing=unit(2, "pt"),
			legend.margin=margin(0, 0, 0, 0),
			axis.text.x=element_text(angle=0, hjust=.5),
			axis.title=element_text(size=10),
			plot.margin=margin(5, 8, 2, 5)
		)
	)
}


#######################
# Command line args
#######################
sldmc_results_output_dir = args[1]
simulation_results_dir = args[2]
tissue_names_file = args[3]
visualization_dir = args[4]
annotation_name_file = args[5]
simulation_oracle_results_dir = args[6]


######################
# Simulations Visualization
######################

# Load in simulation data
simulation_correlation_summary_df = load_simulation_summary(simulation_results_dir, metric="correlation", n_sims=50, eqtl_ss="489", n_anno="6")
simulation_calibration_summary_df = load_simulation_summary(simulation_results_dir, metric="calibration_slope", n_sims=50, eqtl_ss="489", n_anno="6")
sim_corr_oracle_df = load_in_oracle_result(simulation_oracle_results_dir, 3, n_sims=50, n_anno="6")
sim_calibration_oracle_df = load_in_oracle_result(simulation_oracle_results_dir, 4, n_sims=50, n_anno="6")

# Correlation: estimates from each method against the true simulated correlation
simulation_correlation_plot = make_simulation_estimate_bar_plot(simulation_correlation_summary_df, sim_corr_oracle_df, "Correlation of causal\nand predicted eQTL effects", "#3F88C5")
simulation_correlation_output_file = paste0(visualization_dir, "simulation_correlation_estimates.pdf")
ggsave(simulation_correlation_output_file, simulation_correlation_plot, width=7.2, height=3.6)

# Calibration: estimates from each method against the true simulated regression slope
simulation_calibration_plot = make_simulation_estimate_bar_plot(simulation_calibration_summary_df, sim_calibration_oracle_df, "Calibration of causal\non predicted eQTL effects", "#B85C38", y_accuracy=.1)
simulation_calibration_output_file = paste0(visualization_dir, "simulation_calibration_estimates.pdf")
ggsave(simulation_calibration_output_file, simulation_calibration_plot, width=7.2, height=3.6)

# Joint version
joint_simulation <- plot_grid(simulation_calibration_plot, simulation_correlation_plot, ncol=2, labels=c("a", "b"))
simulation_output_file = paste0(visualization_dir, "simulation_joint_stimates.pdf")
ggsave(simulation_output_file, joint_simulation, width=7.2, height=2.5)


######################
# Real data visualization
######################
# Load in tissue info df
tissue_info_df = read.table(tissue_names_file, header=TRUE, sep="\t")

# Map each annotation to its source (sldsc / gene_set / custom) so the forest plots can be split by
# source. The annotation list is (anno_name, source); the diff-stats file itself has no source column.
annotation_source_df = read.table(annotation_name_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
anno_to_source = setNames(as.character(annotation_source_df[[2]]), as.character(annotation_source_df[[1]]))


# Load in meta-analyzed, default S-LDMC results file
per_tissue_sldmc_default_magnitude_df = per_tissue_sldmc_magnitude_results(sldmc_results_output_dir, tissue_info_df)

# Load in meta-analyzed, default S-LDMC results file
sldmc_default_df = read.table(paste0(sldmc_results_output_dir, "sldmc_results_cross_tissue_meta_analyzed_default_bootstrap_stats.txt"), header=TRUE, sep="\t")

# Load in meta-analyzed, default S-LDMC difference-from-intercept results file
sldmc_default_diff_df = read.table(paste0(sldmc_results_output_dir, "sldmc_results_cross_tissue_meta_analyzed_default_intercept_diff_stats.txt"), header=TRUE, sep="\t")

# Load in meta-analyzed, magnitude-stratified S-LDMC results file
sldmc_magnitude_stratified_df = read.table(paste0(sldmc_results_output_dir, "sldmc_results_cross_tissue_meta_analyzed_magnitude_stratified_bootstrap_stats.txt"), header=TRUE, sep="\t")

########################
# Borzoi magnitude bin correlation and calibration plots showing distribution of per-tissue estimates
#######################
per_tissue_magnitude_correlation_plot = make_per_tissue_magnitude_dot_plot(per_tissue_sldmc_default_magnitude_df, "correlation", "S-LDMC\nCorrelation", "#3F88C5")
per_tissue_magnitude_correlation_output_file = paste0(visualization_dir, "sldmc_per_tissue_magnitude_correlation.pdf")
ggsave(per_tissue_magnitude_correlation_output_file, per_tissue_magnitude_correlation_plot, width=7.2, height=3.3)

per_tissue_magnitude_calibration_plot = make_per_tissue_magnitude_dot_plot(per_tissue_sldmc_default_magnitude_df, "calibration_slope", "S-LDMC\nCalibration coefficient", "#B85C38")
per_tissue_magnitude_calibration_output_file = paste0(visualization_dir, "sldmc_per_tissue_magnitude_calibration.pdf")
ggsave(per_tissue_magnitude_calibration_output_file, per_tissue_magnitude_calibration_plot, width=7.2, height=3.3)


########################
# Borzoi magnitude bin correlation and calibration plots showing distribution of per-tissue estimates for only 5 tissues
# x-axis is magnitude bins, and y-axis is correlation or calibration slope bar plot with seperate color for each tissue (and 95% CI error bars)
# Heart left ventricle: #4C72B0
# Brain Cortex #DD8452
# Liver #55A868
# Whole Blood #C44E52
# Muscle Skeletal #8172B3
#######################
five_tissue_colors = c(
	"Heart_Left_Ventricle"="#4C72B0",
	"Brain_Cortex"="#DD8452",
	"Liver"="#55A868",
	"Whole_Blood"="#C44E52",
	"Muscle_Skeletal"="#8172B3"
)

five_tissue_magnitude_correlation_plot = make_five_tissue_magnitude_bar_plot(per_tissue_sldmc_default_magnitude_df, "correlation", "correlation_se", "S-LDMC\nCorrelation", five_tissue_colors, y_accuracy=.1)
five_tissue_magnitude_correlation_output_file = paste0(visualization_dir, "sldmc_five_tissue_magnitude_correlation.pdf")
ggsave(five_tissue_magnitude_correlation_output_file, five_tissue_magnitude_correlation_plot, width=8.5, height=3.6)

# The lowest magnitude bin's calibration estimates are very noisy, so the y-axis is zoomed to the range
# the other bins occupy (a dotted line marks 1 = perfect calibration) and the clipped CI is noted.
five_tissue_calibration_y_limits = c(-1.2, 1.5)
five_tissue_magnitude_calibration_plot = make_five_tissue_magnitude_bar_plot(per_tissue_sldmc_default_magnitude_df, "calibration_slope", "calibration_slope_se", "S-LDMC\nCalibration coefficient", five_tissue_colors, y_limits=five_tissue_calibration_y_limits, y_accuracy=.5, reference_line=1, truncation_note="CIs in the lowest bin extend beyond the axis")
five_tissue_magnitude_calibration_output_file = paste0(visualization_dir, "sldmc_five_tissue_magnitude_calibration.pdf")
ggsave(five_tissue_magnitude_calibration_output_file, five_tissue_magnitude_calibration_plot, width=8.5, height=3.6)

# Make joint calibration + correlation plot for the five tissues (shared x-axis and shared legend)
joint_five_tissue_magnitude_plot = make_stacked_shared_x_plot(five_tissue_magnitude_calibration_plot, five_tissue_magnitude_correlation_plot, c("c", "d"))
ggsave(paste0(visualization_dir, "sldmc_five_tissue_magnitude_correlation_calibration.pdf"), joint_five_tissue_magnitude_plot, width=7.2, height=4.6)


########################
# Annotation forest plots (difference from intercept): separate plots per source, each with its own FDR
########################
# S-LDSC plot drops flanking / allele-frequency annotations; gene-set plot drops MohammadiVg and FinucaneSEG sets
excluded_sldsc_annotation_regex = "flanking|MAF|allele"
excluded_gene_set_annotation_regex = "^MohammadiVg|^FinucaneSEG"

save_source_forest_plot(sldmc_default_diff_df, anno_to_source, "sldsc", excluded_sldsc_annotation_regex, "present", paste0(visualization_dir, "sldmc_sldsc_annotation_forest.pdf"))
save_source_forest_plot(sldmc_default_diff_df, anno_to_source, "gene_set", excluded_gene_set_annotation_regex, "present", paste0(visualization_dir, "sldmc_gene_set_annotation_forest.pdf"))


########################
# Annotation forest plots conditional on borzoi magnitude bin (stratified data): separate plots per source
########################
# Load in meta-analyzed, magnitude-stratified S-LDMC difference-from-intercept results file
sldmc_magnitude_stratified_diff_df = read.table(paste0(sldmc_results_output_dir, "sldmc_results_cross_tissue_meta_analyzed_magnitude_stratified_intercept_diff_stats.txt"), header=TRUE, sep="\t")

for (magnitude_bin_index in 0:4) {
	stratified_present_category = paste0("magnitude_bin", magnitude_bin_index, "Xpresent")
	save_source_forest_plot(sldmc_magnitude_stratified_diff_df, anno_to_source, "sldsc", excluded_sldsc_annotation_regex, stratified_present_category, paste0(visualization_dir, "sldmc_sldsc_annotation_forest_magnitude_bin", magnitude_bin_index, ".pdf"))
	save_source_forest_plot(sldmc_magnitude_stratified_diff_df, anno_to_source, "gene_set", excluded_gene_set_annotation_regex, stratified_present_category, paste0(visualization_dir, "sldmc_gene_set_annotation_forest_magnitude_bin", magnitude_bin_index, ".pdf"))
}


########################
# Borzoi magnitude bin correlation and calibration plots
########################
magnitude_correlation_plot = make_sldmc_magnitude_bar_plot(sldmc_default_df, "correlation", "S-LDMC\nCorrelation", NULL, "#3F88C5")
magnitude_correlation_output_file = paste0(visualization_dir, "sldmc_borzoi_magnitude_correlation.pdf")
ggsave(magnitude_correlation_output_file, magnitude_correlation_plot, width=7.2, height=3.3)

magnitude_calibration_plot = make_sldmc_magnitude_bar_plot(sldmc_default_df, "calibration_slope", "S-LDMC\nCalibration coefficient", NULL, "#B85C38")
magnitude_calibration_output_file = paste0(visualization_dir, "sldmc_borzoi_magnitude_calibration.pdf")
ggsave(magnitude_calibration_output_file, magnitude_calibration_plot, width=7.2, height=3.3)

# joint calibration + correlation plot for the borzoi magnitude bins
joint_magnitude_plot = plot_grid(magnitude_calibration_plot, magnitude_correlation_plot, nrow=1, labels=c("c", "d"))
ggsave(paste0(visualization_dir, "sldmc_borzoi_magnitude_joint.pdf"), joint_magnitude_plot, width=7.2, height=3.3)


########################
# Finer-grained Borzoi magnitude bin correlation and calibration plots
########################
# Bin edges of the 'borzoi_finer_magnitude_bins' annotation (annotate_variant_gene_pairs.py):
# 0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, Inf
finer_magnitude_bin_labels = c("<0.001", "0.001-0.005", "0.005-0.01", "0.01-0.025", "0.025-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.4", ">=0.4")

finer_magnitude_correlation_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "borzoi_finer_magnitude_bins", "correlation", "S-LDMC\nCorrelation", NULL, "Borzoi magnitude bin", "#3F88C5", finer_magnitude_bin_labels, y_accuracy=.1)
finer_magnitude_correlation_output_file = paste0(visualization_dir, "sldmc_borzoi_finer_magnitude_correlation.pdf")
ggsave(finer_magnitude_correlation_output_file, finer_magnitude_correlation_plot, width=7.2, height=3.3)

finer_magnitude_calibration_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "borzoi_finer_magnitude_bins", "calibration_slope", "S-LDMC\nCalibration", NULL, "Borzoi magnitude bin", "#B85C38", finer_magnitude_bin_labels, y_accuracy=.5)
finer_magnitude_calibration_output_file = paste0(visualization_dir, "sldmc_borzoi_finer_magnitude_calibration.pdf")
ggsave(finer_magnitude_calibration_output_file, finer_magnitude_calibration_plot, width=7.2, height=3.3)

# Joint calibration + correlation plot for the finer magnitude bins (x-axis drawn once, under the
# bottom panel; the bin labels stay angled because they are too long to fit horizontally)
joint_finer_magnitude_plot = make_stacked_shared_x_plot(finer_magnitude_calibration_plot, finer_magnitude_correlation_plot, c("c", "d"), shared_legend=FALSE, horizontal_x_labels=FALSE)
joint_finer_magnitude_output_file = paste0(visualization_dir, "sldmc_borzoi_finer_magnitude_joint.pdf")
ggsave(joint_finer_magnitude_output_file, joint_finer_magnitude_plot, width=7.2, height=4.4)



########################
# Allele frequency bin correlation and calibration plots
########################
af_bin_labels = c("<0.05", "0.05-0.1", "0.1-0.2", "0.2-0.3", ">=0.3")

af_correlation_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "af_bins", "correlation", "S-LDMC\nCorrelation", NULL, "Minor allele frequency bin", "#3F88C5", af_bin_labels)
af_correlation_output_file = paste0(visualization_dir, "sldmc_af_correlation.pdf")
ggsave(af_correlation_output_file, af_correlation_plot, width=7.2, height=3.3)

af_calibration_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "af_bins", "calibration_slope", "S-LDMC\nCalibration coefficient", NULL, "Minor allele frequency bin", "#B85C38", af_bin_labels)
af_calibration_output_file = paste0(visualization_dir, "sldmc_af_calibration.pdf")
ggsave(af_calibration_output_file, af_calibration_plot, width=7.2, height=3.3)


########################
# Distance to TSS bin correlation and calibration plots
########################
dist_to_tss_bin_labels = c("0-1 kb", "1-5 kb", "5-25 kb", "25-50 kb", "50-100 kb")

dist_to_tss_correlation_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "dist_to_tss_bins", "correlation", "S-LDMC\nCorrelation", NULL, "Distance to TSS bin", "#3F88C5", dist_to_tss_bin_labels)
dist_to_tss_correlation_output_file = paste0(visualization_dir, "sldmc_dist_to_tss_correlation.pdf")
ggsave(dist_to_tss_correlation_output_file, dist_to_tss_correlation_plot, width=7.2, height=3.3)

dist_to_tss_calibration_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "dist_to_tss_bins", "calibration_slope", "S-LDMC\nCalibration coefficient", NULL, "Distance to TSS bin", "#B85C38", dist_to_tss_bin_labels)
dist_to_tss_calibration_output_file = paste0(visualization_dir, "sldmc_dist_to_tss_calibration.pdf")
ggsave(dist_to_tss_calibration_output_file, dist_to_tss_calibration_plot, width=7.2, height=3.3)


########################
# Distance to TSS bin correlation and calibration plots, stratified by borzoi magnitude bin (facet rows)
########################
dist_to_tss_stratified_correlation_plot = make_sldmc_magnitude_stratified_bar_plot(sldmc_magnitude_stratified_df, "dist_to_tss_bins", "correlation", "S-LDMC\nCorrelation", "Distance to TSS bin", "#3F88C5", dist_to_tss_bin_labels)
dist_to_tss_stratified_correlation_output_file = paste0(visualization_dir, "sldmc_dist_to_tss_magnitude_stratified_correlation.pdf")
ggsave(dist_to_tss_stratified_correlation_output_file, dist_to_tss_stratified_correlation_plot, width=8.4, height=3.6)

dist_to_tss_stratified_calibration_plot = make_sldmc_magnitude_stratified_bar_plot(sldmc_magnitude_stratified_df, "dist_to_tss_bins", "calibration_slope", "S-LDMC\nCalibration coefficient", "Distance to TSS bin", "#B85C38", dist_to_tss_bin_labels)
dist_to_tss_stratified_calibration_output_file = paste0(visualization_dir, "sldmc_dist_to_tss_magnitude_stratified_calibration.pdf")
ggsave(dist_to_tss_stratified_calibration_output_file, dist_to_tss_stratified_calibration_plot, width=8.4, height=3.6)


########################
# Borzoi effect size bin correlation and calibration plots
########################
effect_size_bin_labels = c("< -0.2", "-0.2 to -0.075", "-0.075 to -0.01", "-0.01 to -0.001", "-0.001 to 0", "0 to 0.001", "0.001 to 0.01", "0.01 to 0.075", "0.075 to 0.2", ">= 0.2")

effect_size_correlation_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "borzoi_effect_size_bins", "correlation", "S-LDMC\nCorrelation", NULL, "Borzoi effect size bin", "#3F88C5", effect_size_bin_labels, y_accuracy=.1)
effect_size_correlation_output_file = paste0(visualization_dir, "sldmc_borzoi_effect_size_correlation.pdf")
ggsave(effect_size_correlation_output_file, effect_size_correlation_plot, width=7.2, height=3.3)

effect_size_calibration_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "borzoi_effect_size_bins", "calibration_slope", "S-LDMC\nCalibration", NULL, "Borzoi effect size bin", "#B85C38", effect_size_bin_labels, y_accuracy=.5)
effect_size_calibration_output_file = paste0(visualization_dir, "sldmc_borzoi_effect_size_calibration.pdf")
ggsave(effect_size_calibration_output_file, effect_size_calibration_plot, width=7.2, height=3.3)


# Joint calibration + correlation plot for the borzoi effect size bins (x-axis drawn once, under the
# bottom panel; the bin labels stay angled because they are too long to fit horizontally)
joint_effect_size_plot = make_stacked_shared_x_plot(effect_size_calibration_plot, effect_size_correlation_plot, c("c", "d"), shared_legend=FALSE, horizontal_x_labels=FALSE)
joint_effect_size_output_file = paste0(visualization_dir, "sldmc_borzoi_effect_size_joint.pdf")
ggsave(joint_effect_size_output_file, joint_effect_size_plot, width=7.2, height=4.4)


########################
# Borzoi finer effect size bin correlation and calibration plots
########################
# Bin edges of the 'borzoi_finer_effect_size_bins' annotation (annotate_variant_gene_pairs.py):
# -Inf, -0.4, -0.2, -0.1, -0.05, -0.025, -0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, Inf
finer_effect_size_bin_labels = c("< -0.4", "-0.4 to -0.2", "-0.2 to -0.1", "-0.1 to -0.05", "-0.05 to -0.025", "-0.025 to -0.01", "-0.01 to -0.005", "-0.005 to -0.001", "-0.001 to 0", "0 to 0.001", "0.001 to 0.005", "0.005 to 0.01", "0.01 to 0.025", "0.025 to 0.05", "0.05 to 0.1", "0.1 to 0.2", "0.2 to 0.4", ">= 0.4")

# 18 bins of long labels do not fit in the 7.2in width the coarser plots use, so these are drawn wider
finer_effect_size_correlation_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "borzoi_finer_effect_size_bins", "correlation", "S-LDMC\nCorrelation", NULL, "Borzoi effect size bin", "#3F88C5", finer_effect_size_bin_labels, y_accuracy=.1)
finer_effect_size_correlation_output_file = paste0(visualization_dir, "sldmc_borzoi_finer_effect_size_correlation.pdf")
ggsave(finer_effect_size_correlation_output_file, finer_effect_size_correlation_plot, width=10.5, height=3.6)

finer_effect_size_calibration_plot = make_sldmc_bin_bar_plot(sldmc_default_df, "borzoi_finer_effect_size_bins", "calibration_slope", "S-LDMC\nCalibration coefficient", NULL, "Borzoi effect size bin", "#B85C38", finer_effect_size_bin_labels, y_accuracy=.5)
finer_effect_size_calibration_output_file = paste0(visualization_dir, "sldmc_borzoi_finer_effect_size_calibration.pdf")
ggsave(finer_effect_size_calibration_output_file, finer_effect_size_calibration_plot, width=10.5, height=3.6)

# Joint calibration + correlation plot for the finer borzoi effect size bins (x-axis drawn once, under
# the bottom panel; the bin labels stay angled because they are too long to fit horizontally)
joint_finer_effect_size_plot = make_stacked_shared_x_plot(finer_effect_size_calibration_plot, finer_effect_size_correlation_plot, c("c", "d"), shared_legend=FALSE, horizontal_x_labels=FALSE)
joint_finer_effect_size_output_file = paste0(visualization_dir, "sldmc_borzoi_finer_effect_size_joint.pdf")
ggsave(joint_finer_effect_size_output_file, joint_finer_effect_size_plot, width=10.5, height=5.0)





######################
# Main figure2
#######################

figure2 <- plot_grid(joint_simulation,NULL, joint_effect_size_plot, nrow=3, rel_heights=c(1,.1,1.5))
figure2_output_file = paste0(visualization_dir, "figure2_v1.pdf")
ggsave(figure2_output_file, figure2, width=7.2, height=6.0)

figure2 <- plot_grid(joint_simulation,NULL, joint_finer_magnitude_plot, nrow=3, rel_heights=c(1,.1,1.5))
figure2_output_file = paste0(visualization_dir, "figure2_v2.pdf")
ggsave(figure2_output_file, figure2, width=7.2, height=6.0)