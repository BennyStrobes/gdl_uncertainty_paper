args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(warn=1)

# The inference results directory holds one bootstrap stats file per method per simulation, and the
# two methods do not share a schema: SLDMC reports (annotation_name, category_name) with the
# simulated annotations living in category_name, while the fine-mapped-SNP analysis reports a single
# annotation_name that already matches the annoN names used by the true simulated summaries.
# The directory also holds results from the multi-eQTL-sample-size runs (ss 200/350/450) under the
# same naming scheme, so the patterns pin the sample size to keep those out of these figures.
eqtl_sample_size <- "489"
methods_list <- list(
	list(key="sldmc", label="SLDMC", pattern=paste0("_sim_eqtl_ss_", eqtl_sample_size, "_.*_ld_corr_results_bootstrap_stats\\.txt$")),
	list(key="fine_mapped", label="Fine-mapped SNPs", pattern=paste0("_sim_eqtl_ss_", eqtl_sample_size, "_.*_fm_corr_results_bootstrap_stats\\.txt$"))
)
sldmc_annotation_name <- "sim_signal_prop"

figure_theme <- function() {
	return(
		theme_minimal(base_size=11) +
		theme(
			plot.title=element_text(face="bold", size=12, hjust=.5),
			axis.title=element_text(size=11, color="#25313B"),
			axis.text=element_text(size=10, color="#33424F"),
			panel.grid.minor=element_blank(),
			panel.grid.major.x=element_blank(),
			panel.grid.major.y=element_blank(),
			panel.background=element_rect(fill="#FBFCFD", color=NA),
			plot.background=element_rect(fill="#FBFCFD", color=NA),
			axis.line.x=element_line(colour="#7E8A97", linewidth=.35),
			axis.line.y=element_line(colour="#7E8A97", linewidth=.35),
			legend.position="none"
		)
	)
}

# Annotations arrive as annoN names but are displayed as numbered simulated categories (1, 2, ...),
# following the order of the supplied factor levels.
annotation_x_scale <- function(annotation_levels) {
	return(scale_x_discrete(labels=setNames(as.character(seq_along(annotation_levels)), annotation_levels)))
}
annotation_x_label <- "Simulated annotation category"

make_bias_barplot <- function(summary_df, plot_title) {
	summary_df$annotation_name <- factor(summary_df$annotation_name, levels=summary_df$annotation_name)
	summary_df$xpos <- seq_len(nrow(summary_df))
	return(
		ggplot(summary_df, aes(x=annotation_name, y=mean)) +
			geom_col(fill="#688B84", color="#365854", width=.68, linewidth=.35) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.16, linewidth=.45, color="#22313B") +
			geom_segment(aes(x=xpos-.30, xend=xpos+.30, y=true_simulated_value, yend=true_simulated_value), color="#D06A4B", linetype="dashed", linewidth=.75, inherit.aes=FALSE) +
			geom_point(aes(y=true_simulated_value), color="#D06A4B", size=1.8, inherit.aes=TRUE) +
			annotation_x_scale(levels(summary_df$annotation_name)) +
			xlab(annotation_x_label) +
			ylab("Value") +
			ggtitle(plot_title) +
			figure_theme()
	)
}

make_coverage_plot <- function(coverage_df, plot_title) {
	coverage_df$annotation_name <- factor(coverage_df$annotation_name, levels=unique(coverage_df$annotation_name))
	coverage_df$coverage_label <- factor(coverage_df$coverage_label, levels=c("99%", "95%", "90%"))
	y_min <- min(c(coverage_df$ci_lower, coverage_df$nominal_coverage), na.rm=TRUE)
	y_max <- max(c(coverage_df$ci_upper, coverage_df$nominal_coverage), na.rm=TRUE)
	y_padding <- 0.02
	y_min <- max(0, y_min - y_padding)
	y_max <- min(1, y_max + y_padding)
	dodge_width <- 0.35
	return(
		ggplot(coverage_df, aes(x=annotation_name, y=observed_coverage, color=coverage_label, group=coverage_label)) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.14, linewidth=.45, position=position_dodge(width=dodge_width)) +
			geom_line(aes(group=coverage_label), linewidth=.7, position=position_dodge(width=dodge_width)) +
			geom_point(size=2, position=position_dodge(width=dodge_width)) +
			geom_hline(aes(yintercept=nominal_coverage, color=coverage_label), linetype="dotted", linewidth=.6, show.legend=FALSE) +
			scale_color_manual(values=c("90%"="#4C78A8", "95%"="#D06A4B", "99%"="#7A5EA6"), name="") +
			annotation_x_scale(levels(coverage_df$annotation_name)) +
			xlab(annotation_x_label) +
			ylab("Coverage") +
			ggtitle(plot_title) +
			coord_cartesian(ylim=c(y_min, y_max)) +
			figure_theme() +
			theme(legend.position="right")
	)
}

make_se_comparison_plot <- function(se_df, plot_title) {
	long_df <- rbind(
		data.frame(
			annotation_name=se_df$annotation_name,
			se_type="Average estimated SE",
			value=se_df$mean_estimated_se,
			ci_lower=se_df$mean_estimated_se_ci_lower,
			ci_upper=se_df$mean_estimated_se_ci_upper,
			stringsAsFactors=FALSE
		),
		data.frame(
			annotation_name=se_df$annotation_name,
			se_type="Empirical SE",
			value=se_df$empirical_se,
			ci_lower=se_df$empirical_se_ci_lower,
			ci_upper=se_df$empirical_se_ci_upper,
			stringsAsFactors=FALSE
		)
	)
	long_df$annotation_name <- factor(long_df$annotation_name, levels=se_df$annotation_name)
	long_df$se_type <- factor(long_df$se_type, levels=c("Average estimated SE", "Empirical SE"))
	dodge_width <- 0.7
	return(
		ggplot(long_df, aes(x=annotation_name, y=value, fill=se_type, group=se_type)) +
			geom_col(width=.62, color="#22313B", linewidth=.35, position=position_dodge(width=dodge_width)) +
			geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=.16, linewidth=.45, color="#22313B", position=position_dodge(width=dodge_width)) +
			scale_fill_manual(values=c("Average estimated SE"="#3E7A34", "Empirical SE"="#B863A0"), name="") +
			annotation_x_scale(levels(long_df$annotation_name)) +
			xlab(annotation_x_label) +
			ylab("Standard error") +
			ggtitle(plot_title) +
			figure_theme() +
			theme(legend.position="right")
	)
}

# Ratio of average estimated SE to empirical SE, on a log axis so that a 25% overestimate and a 20%
# underestimate sit the same distance from the calibrated line at 1.
make_se_ratio_plot <- function(se_df, plot_title) {
	se_df$annotation_name <- factor(se_df$annotation_name, levels=se_df$annotation_name)
	ratio_plot <- ggplot(se_df, aes(x=annotation_name, y=se_ratio, group=1)) +
		geom_hline(yintercept=1, linetype="dashed", linewidth=.6, color="#D06A4B") +
		geom_errorbar(aes(ymin=se_ratio_ci_lower, ymax=se_ratio_ci_upper), width=.14, linewidth=.45, color="#33424F") +
		geom_line(linewidth=.7, color="#4C78A8") +
		geom_point(size=2, color="#22313B") +
		annotation_x_scale(levels(se_df$annotation_name)) +
		xlab(annotation_x_label) +
		ylab("Average estimated SE / empirical SE") +
		ggtitle(plot_title) +
		figure_theme()

	ratio_values <- c(se_df$se_ratio, se_df$se_ratio_ci_lower, se_df$se_ratio_ci_upper)
	if (all(is.finite(ratio_values)) && all(ratio_values > 0)) {
		ratio_plot <- ratio_plot + scale_y_log10()
	}
	return(ratio_plot)
}

extract_simulation_iter_from_file <- function(file_name) {
	matches <- regmatches(file_name, regexpr("sim[0-9]+", file_name))
	if (length(matches) == 0 || matches[1] == "") {
		return(NA)
	}
	return(matches[1])
}

order_annotation_names <- function(annotation_names) {
	annotation_names <- unique(annotation_names)
	if (length(annotation_names) > 0 && all(grepl("^anno[0-9]+$", annotation_names))) {
		return(annotation_names[order(as.numeric(sub("^anno", "", annotation_names)))])
	}
	return(annotation_names[order(annotation_names)])
}

# Read every bootstrap stats file matching a single method's naming convention, keeping only the
# columns shared across schemas so that files from different methods are never pooled together.
read_bootstrap_stats_files <- function(results_dir, file_pattern) {
	result_files <- list.files(results_dir, pattern=file_pattern, full.names=TRUE)
	if (length(result_files) == 0) {
		return(NULL)
	}

	kept_columns <- c("simulation_iter", "annotation_name", "category_name", "output_name", "mean", "bootstrap_se")
	all_est <- data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp <- read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		required_columns <- c("annotation_name", "output_name", "mean", "bootstrap_se")
		missing_columns <- required_columns[!(required_columns %in% colnames(tmp))]
		if (length(missing_columns) > 0) {
			stop(paste("Missing column(s)", paste(missing_columns, collapse=", "), "in", result_file))
		}
		if (!("category_name" %in% colnames(tmp))) {
			tmp$category_name <- NA
		}
		tmp$simulation_iter <- extract_simulation_iter_from_file(basename(result_file))
		all_est <- rbind(all_est, tmp[, kept_columns])
	}
	if (nrow(all_est) == 0) {
		return(NULL)
	}
	return(all_est)
}

# SLDMC names its output categories 'prop_<gaussian proportion>', so recovering the annoN names used
# everywhere else means going through the category file, where category_index is the one-hot index of
# the corresponding annoN column.
load_sldmc_category_map <- function(correlation_borzoi_est_effect_dir) {
	category_files <- list.files(correlation_borzoi_est_effect_dir, pattern="_sldmc_categories\\.txt$", full.names=TRUE)
	category_map <- data.frame()
	for (category_file in category_files) {
		if (is.na(file.info(category_file)$size) || file.info(category_file)$size == 0) {
			next
		}
		tmp <- read.table(category_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp <- tmp[tmp$anno_name == sldmc_annotation_name, ]
		if (nrow(tmp) == 0) {
			next
		}
		category_map <- rbind(
			category_map,
			data.frame(
				category_name=tmp$category_name,
				mapped_annotation_name=paste0("anno", tmp$category_index),
				stringsAsFactors=FALSE
			)
		)
	}
	if (nrow(category_map) == 0) {
		return(NULL)
	}
	return(unique(category_map))
}

# Fallback when no category file is available: the simulated proportions are generated on an
# ascending grid, so sorting the prop_ categories numerically recovers the annotation index.
build_sldmc_category_map_from_names <- function(category_names) {
	category_names <- unique(category_names)
	if (!all(grepl("^prop_", category_names))) {
		return(NULL)
	}
	prop_values <- suppressWarnings(as.numeric(sub("^prop_", "", category_names)))
	if (any(is.na(prop_values))) {
		return(NULL)
	}
	category_names <- category_names[order(prop_values)]
	return(data.frame(
		category_name=category_names,
		mapped_annotation_name=paste0("anno", seq_along(category_names) - 1),
		stringsAsFactors=FALSE
	))
}

harmonize_sldmc_annotation_names <- function(all_est, correlation_borzoi_est_effect_dir) {
	# Result files predating the SLDMC refactor have no category_name and so cannot be placed on the
	# same annotation axis; drop them rather than pooling two different estimators.
	n_legacy_rows <- sum(is.na(all_est$category_name))
	if (n_legacy_rows > 0) {
		message(paste("Ignoring", n_legacy_rows, "SLDMC result row(s) with no category_name (pre-SLDMC output format)"))
	}
	all_est <- all_est[!is.na(all_est$category_name) & all_est$annotation_name == sldmc_annotation_name, ]
	if (nrow(all_est) == 0) {
		message(paste("No SLDMC result rows with annotation_name", sldmc_annotation_name, "were found"))
		return(NULL)
	}
	category_map <- load_sldmc_category_map(correlation_borzoi_est_effect_dir)
	if (is.null(category_map)) {
		category_map <- build_sldmc_category_map_from_names(all_est$category_name)
	}
	if (is.null(category_map)) {
		stop(paste("Could not map SLDMC category names onto annotation names using category files in", correlation_borzoi_est_effect_dir))
	}
	merged_df <- merge(all_est, category_map, by="category_name", all.x=TRUE)
	unmapped_categories <- unique(merged_df$category_name[is.na(merged_df$mapped_annotation_name)])
	if (length(unmapped_categories) > 0) {
		stop(paste("Unmapped SLDMC category name(s):", paste(unmapped_categories, collapse=", ")))
	}
	merged_df$annotation_name <- merged_df$mapped_annotation_name
	merged_df$mapped_annotation_name <- NULL
	return(merged_df)
}

load_estimated_corr_per_sim <- function(correlation_inference_results_dir, correlation_borzoi_est_effect_dir, method) {
	all_est <- read_bootstrap_stats_files(correlation_inference_results_dir, method$pattern)
	if (is.null(all_est)) {
		return(NULL)
	}
	if (method$key == "sldmc") {
		all_est <- harmonize_sldmc_annotation_names(all_est, correlation_borzoi_est_effect_dir)
		if (is.null(all_est)) {
			return(NULL)
		}
	}
	duplicate_keys <- paste(all_est$simulation_iter, all_est$annotation_name, all_est$output_name, sep=":")
	if (any(duplicated(duplicate_keys))) {
		stop(paste("Duplicate (simulation, annotation, output) rows found for method", method$key))
	}
	return(all_est)
}

summarize_estimated_corr_across_sims <- function(all_est, output_name) {
	subset_df <- all_est[all_est$output_name == output_name, ]
	annotation_names <- order_annotation_names(subset_df$annotation_name)

	output_summary <- data.frame(
		annotation_name=annotation_names,
		mean=rep(NA, length(annotation_names)),
		ci_lower=rep(NA, length(annotation_names)),
		ci_upper=rep(NA, length(annotation_names)),
		n_sims=rep(NA, length(annotation_names)),
		stringsAsFactors=FALSE
	)

	for (anno_iter in seq_along(annotation_names)) {
		anno_name <- annotation_names[anno_iter]
		anno_df <- subset_df[subset_df$annotation_name == anno_name, ]
		anno_means <- anno_df$mean
		n_sims <- length(anno_means)
		mean_val <- mean(anno_means)
		sem_val <- sd(anno_means)/sqrt(n_sims)

		output_summary$mean[anno_iter] <- mean_val
		output_summary$ci_lower[anno_iter] <- mean_val - 1.96*sem_val
		output_summary$ci_upper[anno_iter] <- mean_val + 1.96*sem_val
		output_summary$n_sims[anno_iter] <- n_sims
	}
	return(output_summary)
}

load_true_corr_per_sim <- function(correlation_borzoi_est_effect_dir) {
	result_files <- list.files(correlation_borzoi_est_effect_dir, pattern="_true_sim_effect_summary\\.txt$", full.names=TRUE)
	if (length(result_files) == 0) {
		stop(paste("No true simulated summary files found in", correlation_borzoi_est_effect_dir))
	}

	all_true <- data.frame()
	for (result_file in result_files) {
		if (is.na(file.info(result_file)$size) || file.info(result_file)$size == 0) {
			next
		}
		tmp <- read.table(result_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		tmp$simulation_iter <- extract_simulation_iter_from_file(basename(result_file))
		all_true <- rbind(all_true, tmp)
	}
	if (nrow(all_true) == 0) {
		stop(paste("True simulated summary files were found in", correlation_borzoi_est_effect_dir, "but all were empty"))
	}
	return(all_true)
}

summarize_true_corr_across_sims <- function(all_true, true_column_name) {
	annotation_names <- order_annotation_names(all_true$annotation_name)
	true_df <- data.frame(
		annotation_name=annotation_names,
		true_simulated_value=rep(NA, length(annotation_names)),
		stringsAsFactors=FALSE
	)
	for (anno_iter in seq_along(annotation_names)) {
		anno_df <- all_true[all_true$annotation_name == annotation_names[anno_iter], ]
		true_df$true_simulated_value[anno_iter] <- mean(anno_df[, true_column_name])
	}
	return(true_df)
}

merge_estimated_and_true_summary <- function(estimated_summary_df, true_summary_df) {
	merged_df <- merge(estimated_summary_df, true_summary_df, by="annotation_name", all.x=TRUE, sort=FALSE)
	merged_df <- merged_df[match(estimated_summary_df$annotation_name, merged_df$annotation_name), ]
	rownames(merged_df) <- NULL
	return(merged_df)
}

make_coverage_summary <- function(all_est, all_true, output_name, true_column_name) {
	est_df <- all_est[all_est$output_name == output_name, c("simulation_iter", "annotation_name", "mean", "bootstrap_se")]
	true_df <- all_true[, c("simulation_iter", "annotation_name", true_column_name)]
	colnames(true_df)[3] <- "true_value"
	merged_df <- merge(est_df, true_df, by=c("simulation_iter", "annotation_name"))
	if (nrow(merged_df) == 0) {
		stop(paste("No overlapping (simulation, annotation) pairs between estimated and true results for", output_name))
	}
	coverage_levels <- c(.90, .95, .99)
	coverage_summary <- data.frame()
	for (coverage_level in coverage_levels) {
		z_value <- qnorm((1.0 + coverage_level)/2.0)
		tmp_df <- merged_df
		tmp_df$ci_lower <- tmp_df$mean - z_value*tmp_df$bootstrap_se
		tmp_df$ci_upper <- tmp_df$mean + z_value*tmp_df$bootstrap_se
		tmp_df$covered <- (tmp_df$true_value >= tmp_df$ci_lower) & (tmp_df$true_value <= tmp_df$ci_upper)
		annotation_names <- order_annotation_names(tmp_df$annotation_name)
		for (anno_name in annotation_names) {
			anno_df <- tmp_df[tmp_df$annotation_name == anno_name, ]
			observed_coverage <- mean(anno_df$covered)
			n_sims <- nrow(anno_df)
			coverage_sem <- sqrt(observed_coverage*(1.0-observed_coverage)/n_sims)
			coverage_ci_lower <- max(0, observed_coverage - 1.96*coverage_sem)
			coverage_ci_upper <- min(1, observed_coverage + 1.96*coverage_sem)
			coverage_summary <- rbind(
				coverage_summary,
				data.frame(
					annotation_name=anno_name,
					coverage_label=paste0(as.integer(coverage_level*100), "%"),
					nominal_coverage=coverage_level,
					observed_coverage=observed_coverage,
					coverage_sem=coverage_sem,
					ci_lower=coverage_ci_lower,
					ci_upper=coverage_ci_upper,
					stringsAsFactors=FALSE
				)
			)
		}
	}
	return(coverage_summary)
}

# Empirical SE is sd() of the residuals rather than their root-mean-square, so it measures the spread
# of the estimator around its own mean and is not inflated by bias. Estimated SEs are averaged in
# variance space (root-mean-square) so the comparison to that spread is variance against variance.
compute_se_calibration_stats <- function(residuals, estimated_ses) {
	empirical_se <- sd(residuals)
	mean_estimated_se <- sqrt(mean(estimated_ses^2))
	if (!is.finite(empirical_se) || empirical_se == 0) {
		se_ratio <- NA
	} else {
		se_ratio <- mean_estimated_se/empirical_se
	}
	return(c(empirical_se=empirical_se, mean_estimated_se=mean_estimated_se, se_ratio=se_ratio))
}

# Resampling whole simulations keeps each residual paired with the estimated SE from the same
# simulation, so the ratio's interval accounts for both quantities moving together.
bootstrap_se_calibration_stats <- function(residuals, estimated_ses, n_bootstraps) {
	n_sims <- length(residuals)
	bootstrap_draws <- matrix(NA, nrow=n_bootstraps, ncol=3)
	colnames(bootstrap_draws) <- c("empirical_se", "mean_estimated_se", "se_ratio")
	for (bootstrap_iter in seq_len(n_bootstraps)) {
		bootstrap_indices <- sample(seq_len(n_sims), size=n_sims, replace=TRUE)
		bootstrap_draws[bootstrap_iter, ] <- compute_se_calibration_stats(residuals[bootstrap_indices], estimated_ses[bootstrap_indices])
	}
	return(bootstrap_draws)
}

bootstrap_ci <- function(bootstrap_values) {
	if (all(is.na(bootstrap_values))) {
		return(c(NA, NA))
	}
	return(unname(quantile(bootstrap_values, c(.025, .975), na.rm=TRUE)))
}

make_se_calibration_summary <- function(all_est, all_true, output_name, true_column_name, n_bootstraps=1000) {
	est_df <- all_est[all_est$output_name == output_name, c("simulation_iter", "annotation_name", "mean", "bootstrap_se")]
	true_df <- all_true[, c("simulation_iter", "annotation_name", true_column_name)]
	colnames(true_df)[3] <- "true_value"
	merged_df <- merge(est_df, true_df, by=c("simulation_iter", "annotation_name"))

	n_dropped_rows <- sum(!(is.finite(merged_df$mean) & is.finite(merged_df$bootstrap_se) & is.finite(merged_df$true_value)))
	if (n_dropped_rows > 0) {
		message(paste("Dropping", n_dropped_rows, "non-finite (simulation, annotation) row(s) for", output_name))
	}
	merged_df <- merged_df[is.finite(merged_df$mean) & is.finite(merged_df$bootstrap_se) & is.finite(merged_df$true_value), ]
	if (nrow(merged_df) == 0) {
		stop(paste("No usable (simulation, annotation) pairs for standard error calibration of", output_name))
	}
	merged_df$residual <- merged_df$mean - merged_df$true_value

	annotation_names <- order_annotation_names(merged_df$annotation_name)
	se_summary <- data.frame()
	for (anno_name in annotation_names) {
		anno_df <- merged_df[merged_df$annotation_name == anno_name, ]
		n_sims <- nrow(anno_df)
		if (n_sims < 2) {
			message(paste("Annotation", anno_name, "has fewer than 2 usable simulations for", output_name, "- skipping"))
			next
		}
		observed_stats <- compute_se_calibration_stats(anno_df$residual, anno_df$bootstrap_se)
		bootstrap_draws <- bootstrap_se_calibration_stats(anno_df$residual, anno_df$bootstrap_se, n_bootstraps)
		empirical_ci <- bootstrap_ci(bootstrap_draws[, "empirical_se"])
		estimated_ci <- bootstrap_ci(bootstrap_draws[, "mean_estimated_se"])
		ratio_ci <- bootstrap_ci(bootstrap_draws[, "se_ratio"])

		se_summary <- rbind(
			se_summary,
			data.frame(
				annotation_name=anno_name,
				n_sims=n_sims,
				mean_residual=mean(anno_df$residual),
				empirical_se=observed_stats["empirical_se"],
				empirical_se_ci_lower=empirical_ci[1],
				empirical_se_ci_upper=empirical_ci[2],
				mean_estimated_se=observed_stats["mean_estimated_se"],
				mean_estimated_se_ci_lower=estimated_ci[1],
				mean_estimated_se_ci_upper=estimated_ci[2],
				se_ratio=observed_stats["se_ratio"],
				se_ratio_ci_lower=ratio_ci[1],
				se_ratio_ci_upper=ratio_ci[2],
				stringsAsFactors=FALSE
			)
		)
	}
	if (nrow(se_summary) == 0) {
		return(NULL)
	}
	rownames(se_summary) <- NULL
	return(se_summary)
}


####################
# Command line args
####################
correlation_inference_results_dir <- args[1]
correlation_borzoi_est_effect_dir <- args[2]
viz_dir <- args[3]

if (!dir.exists(viz_dir)) {
	dir.create(viz_dir, recursive=TRUE)
}

# The standard error calibration intervals are bootstrapped, so fix the seed to keep figures stable
# across reruns on the same results.
set.seed(1)

#################
# Load in true simulated values
################
all_true_results <- load_true_corr_per_sim(correlation_borzoi_est_effect_dir)
true_result_summary <- list(
	correlation=summarize_true_corr_across_sims(all_true_results, "pearson_correlation"),
	calibration_slope=summarize_true_corr_across_sims(all_true_results, "regression_slope")
)


#################
# One set of bias and coverage plots per inference method
################
output_specs <- list(
	list(key="correlation", output_name="correlation", true_column_name="pearson_correlation", label="Correlation"),
	list(key="calibration_slope", output_name="calibration_slope", true_column_name="regression_slope", label="Calibration Slope")
)

n_methods_plotted <- 0
for (method in methods_list) {
	all_estimated_results <- load_estimated_corr_per_sim(correlation_inference_results_dir, correlation_borzoi_est_effect_dir, method)
	if (is.null(all_estimated_results)) {
		message(paste("No usable result files for method", method$key, "in", correlation_inference_results_dir, "- skipping"))
		next
	}
	n_methods_plotted <- n_methods_plotted + 1
	se_plots_by_output <- list()

	for (output_spec in output_specs) {
		if (!(output_spec$output_name %in% all_estimated_results$output_name)) {
			message(paste("Method", method$key, "has no", output_spec$output_name, "rows - skipping"))
			next
		}

		#################
		# Bias plot
		#################
		estimated_summary_df <- summarize_estimated_corr_across_sims(all_estimated_results, output_spec$output_name)
		estimated_summary_df <- merge_estimated_and_true_summary(estimated_summary_df, true_result_summary[[output_spec$key]])
		summary_output_file <- file.path(viz_dir, paste0("estimated_", output_spec$key, "_summary_across_sims_", method$key, ".txt"))
		write.table(estimated_summary_df, file=summary_output_file, sep="\t", quote=FALSE, row.names=FALSE)

		bias_plot <- make_bias_barplot(estimated_summary_df, paste0(output_spec$label, " (", method$label, ")"))
		bias_plot_file <- file.path(viz_dir, paste0("estimated_", output_spec$key, "_bias_barplot_", method$key, ".pdf"))
		ggsave(filename=bias_plot_file, plot=bias_plot, width=7.2, height=3.5)

		#################
		# Coverage plot
		#################
		coverage_df <- make_coverage_summary(all_estimated_results, all_true_results, output_spec$output_name, output_spec$true_column_name)
		coverage_output_file <- file.path(viz_dir, paste0("estimated_", output_spec$key, "_coverage_summary_", method$key, ".txt"))
		write.table(coverage_df, file=coverage_output_file, sep="\t", quote=FALSE, row.names=FALSE)

		coverage_plot <- make_coverage_plot(coverage_df, paste0(output_spec$label, " Coverage (", method$label, ")"))
		coverage_plot_file <- file.path(viz_dir, paste0("estimated_", output_spec$key, "_coverage_plot_", method$key, ".pdf"))
		ggsave(filename=coverage_plot_file, plot=coverage_plot, width=7.2, height=3.6)


		#################
		# Empirical standard error plot
		#################
		se_calibration_df <- make_se_calibration_summary(all_estimated_results, all_true_results, output_spec$output_name, output_spec$true_column_name)
		if (is.null(se_calibration_df)) {
			message(paste("No annotation had enough simulations for standard error calibration of", output_spec$output_name, "-", method$key))
			next
		}
		se_output_file <- file.path(viz_dir, paste0("estimated_", output_spec$key, "_se_calibration_summary_", method$key, ".txt"))
		write.table(se_calibration_df, file=se_output_file, sep="\t", quote=FALSE, row.names=FALSE)

		se_comparison_plot <- make_se_comparison_plot(se_calibration_df, paste0(output_spec$label, " Standard Errors (", method$label, ")"))
		se_comparison_plot_file <- file.path(viz_dir, paste0("estimated_", output_spec$key, "_se_comparison_plot_", method$key, ".pdf"))
		ggsave(filename=se_comparison_plot_file, plot=se_comparison_plot, width=7.2, height=3.6)
		se_plots_by_output[[output_spec$key]] <- se_comparison_plot

	}

	#################
	# Combined standard error plot (calibration slope on top, correlation on bottom)
	#################
	if (!is.null(se_plots_by_output$calibration_slope) && !is.null(se_plots_by_output$correlation)) {
		combined_se_plot <- plot_grid(
			se_plots_by_output$calibration_slope,
			se_plots_by_output$correlation,
			ncol=1, align="v", axis="lr"
		)
		combined_se_plot_file <- file.path(viz_dir, paste0("estimated_combined_se_comparison_plot_", method$key, ".pdf"))
		ggsave(filename=combined_se_plot_file, plot=combined_se_plot, width=7.2, height=6.0)
	}
}

if (n_methods_plotted == 0) {
	stop(paste("No estimated result files found in", correlation_inference_results_dir))
}
