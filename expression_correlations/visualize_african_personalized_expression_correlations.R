args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
options(warn=1)

# European vs African-ancestry comparison of personalized expression prediction.
# Inputs are the per-gene outputs of personalized_expression_correlations_per_tissue.py (European) and
# african_personalized_expression_correlations_per_tissue.py (African ancestry). Both use the same Borzoi
# effects and the same European-trained S-LDMC rescaling; they differ in the evaluation samples (LD,
# allele frequencies, expression). Genes are matched by (tissue, gene_id) and only genes with finite
# predicted and observed correlations in both ancestries are compared.
# (The single-ancestry AFR plots are made by running visualize_personalized_expression_correlations.R on
# the AFR results directory.)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


# Heritable genes: European cis-SNP heritability LRT p-value below this threshold. Selected on the European
# results only, so the African observed correlations are not subject to winner's curse from the selection.
heritability_pvalue_threshold = 0.05

# Observed correlation between the rescaled prediction and observed expression
observed_r_col = "rescaled_expression_correlation"

population_levels = c("EUR", "AFR")
quantity_colors = c("Predicted r"="#6B7280", "Observed r"="#B91C1C")
population_colors = c("EUR"="#2563EB", "AFR"="#D97706")
population_legend_labels = c("EUR"="European", "AFR"="African")
reference_color = "#3F3F46"

# Expression-FSR bins (as in visualize_personalized_expression_correlations.R)
fsr_bin_breaks = c(0, 0.1, 0.2, 0.3, 0.4, Inf)
fsr_bin_labels = c("[0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", ">0.4")

# Predefined predicted-correlation bins (as in visualize_personalized_expression_correlations.R)
predicted_r_bin_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.6, Inf)
predicted_r_bin_labels = c("[0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.6]", ">0.6")

# Observed cis-SNP heritability used for observed r / sqrt(h2): the unbiased Haseman-Elston estimate
# (as in visualize_personalized_expression_correlations.R), taken from each ancestry's own results
observed_h2_col = "cis_snp_h2_he"

# Upper limit of the y-axis of the two-row predicted r calibration plot: a small-sample AFR bin whose mean
# observed h2 is close to 0 can give a very large observed r / sqrt(h2), which would otherwise squash
# every other panel (points / CIs beyond the limit are clipped, not dropped; see the bin summary file)
two_row_predicted_r_y_max_cap = 1.5

# Smallest N shown in the top-N-genes plot (as in visualize_personalized_expression_correlations.R)
top_n_min_genes = 10

# Residual-variance models for the predicted correlation with genetic expression in the top-N-genes plot
# ("" = constant within borzoi magnitude bin; "_af_specific" = allele-frequency-specific)
ranking_model_suffixes = c("", "_af_specific")


read_per_tissue_expression_correlation_file <- function(results_file) {
	# Read one per-gene output file, skipping (with a warning) any line that does not have the header's
	# number of columns and a last line with no newline terminator: a job that is still running or died
	# mid-write leaves a cut-off last line. (Same as in visualize_personalized_expression_correlations.R)
	file_text = readChar(results_file, nchars=file.size(results_file), useBytes=TRUE)
	file_lines = strsplit(file_text, "\n", fixed=TRUE)[[1]]
	n_fields_per_line = nchar(gsub("[^\t]", "", file_lines)) + 1
	bad_lines = which(n_fields_per_line != n_fields_per_line[1])
	if (endsWith(file_text, "\n") == FALSE) {
		bad_lines = union(bad_lines, length(file_lines))
	}
	if (length(bad_lines) > 0) {
		print(paste0("WARNING: ", results_file, ": skipping ", length(bad_lines), " cut-off line(s) out of ", length(file_lines), " (line ", paste(head(bad_lines, 5), collapse=", "), "); the job is probably still running or died mid-write"))
		file_lines = file_lines[-bad_lines]
	}
	return(read.table(text=file_lines, header=TRUE, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE))
}


load_per_tissue_expression_correlations <- function(per_tissue_expression_correlation_dir, tissue_info_df) {
	# Stack every tissue's per-gene output file
	# (<tissue>_<sample>_personalized_expression_prediction.txt) into one long data frame with a
	# target_tissue column. Tissues without a results file are skipped with a warning.
	# (Same as in visualize_personalized_expression_correlations.R)
	results_df = data.frame()
	for (row_iter in seq_len(nrow(tissue_info_df))) {
		target_tissue = tissue_info_df$GTEx_tissue[row_iter]
		target_sample = tissue_info_df$target_identifier[row_iter]
		results_file = paste0(per_tissue_expression_correlation_dir, target_tissue, "_", target_sample, "_personalized_expression_prediction.txt")
		if (file.exists(results_file) == FALSE) {
			print(paste("WARNING: per-tissue expression correlation file not found; skipping:", results_file))
			next
		}
		tissue_df = read_per_tissue_expression_correlation_file(results_file)
		tissue_df$target_tissue = target_tissue
		results_df = rbind(results_df, tissue_df)
	}
	if (nrow(results_df) == 0) {
		stop(paste("No per-tissue expression correlation files found in", per_tissue_expression_correlation_dir))
	}
	return(results_df)
}


merge_eur_afr_results <- function(eur_df, afr_df) {
	# One row per (tissue, gene) analyzed in both ancestries, with _eur / _afr columns. Adds the predicted
	# correlation with observed expression, sqrt(Var_i(X mu)), computed in each ancestry's own samples, and
	# keeps only genes with finite predicted and observed correlations in both.
	keep_cols = c("target_tissue", "gene_id", observed_r_col, "rescaled_predicted_expression_variance", "predicted_cis_snp_h2", "predicted_cis_snp_h2_af_specific", "expression_FSR_af_specific", "cis_snp_h2_pvalue", "n_samples", "n_cis_variants")
	merged_df = merge(eur_df[, keep_cols], afr_df[, keep_cols], by=c("target_tissue", "gene_id"), suffixes=c("_eur", "_afr"))
	for (suffix in c("_eur", "_afr")) {
		merged_df[[paste0("predicted_r", suffix)]] = sqrt(merged_df[[paste0("rescaled_predicted_expression_variance", suffix)]])
		merged_df[[paste0("observed_r", suffix)]] = merged_df[[paste0(observed_r_col, suffix)]]
		# Expected correlation with genetic expression (AF-specific): sqrt(Var_i(X mu) / predicted cis-SNP h2 (AF-specific))
		merged_df[[paste0("predicted_genetic_r_af_specific", suffix)]] = sqrt(merged_df[[paste0("rescaled_predicted_expression_variance", suffix)]]/merged_df[[paste0("predicted_cis_snp_h2_af_specific", suffix)]])
	}
	complete = is.finite(merged_df$predicted_r_eur) & is.finite(merged_df$predicted_r_afr) & is.finite(merged_df$observed_r_eur) & is.finite(merged_df$observed_r_afr)
	print(paste("Genes analyzed in both ancestries:", nrow(merged_df), "; with finite predicted and observed correlations in both:", sum(complete)))
	return(merged_df[complete, ])
}


compute_tissue_sample_size_df <- function(eur_df, afr_df, merged_df, tissue_colors) {
	# Per tissue: expression sample sizes and gene counts in each ancestry, genes compared in both, and
	# genes compared in both that are heritable in Europeans
	rows = list()
	for (target_tissue in names(tissue_colors)) {
		eur_tissue_df = eur_df[eur_df$target_tissue == target_tissue, ]
		afr_tissue_df = afr_df[afr_df$target_tissue == target_tissue, ]
		merged_tissue_df = merged_df[merged_df$target_tissue == target_tissue, ]
		rows[[length(rows) + 1]] = data.frame(
			target_tissue=target_tissue,
			n_samples_eur=ifelse(nrow(eur_tissue_df) > 0, max(eur_tissue_df$n_samples), NA),
			n_samples_afr=ifelse(nrow(afr_tissue_df) > 0, max(afr_tissue_df$n_samples), NA),
			n_genes_eur=nrow(eur_tissue_df),
			n_genes_afr=nrow(afr_tissue_df),
			n_genes_both=nrow(merged_tissue_df),
			n_genes_both_heritable_eur=sum(!is.na(merged_tissue_df$cis_snp_h2_pvalue_eur) & merged_tissue_df$cis_snp_h2_pvalue_eur < heritability_pvalue_threshold),
			stringsAsFactors=FALSE
		)
	}
	return(do.call(rbind, rows))
}


make_tissue_labels <- function(sample_size_df, tissue_colors) {
	# Facet labels: tissue name with the EUR and AFR expression sample sizes
	tissue_labels = c()
	for (target_tissue in names(tissue_colors)) {
		row = sample_size_df[sample_size_df$target_tissue == target_tissue, ]
		tissue_labels[target_tissue] = paste0(gsub("_", " ", target_tissue), "\n(EUR n=", row$n_samples_eur, ", AFR n=", row$n_samples_afr, ")")
	}
	return(tissue_labels)
}


compute_population_mean_r_summary_df <- function(merged_df, tissue_colors, tissue_labels) {
	# One row per (tissue, ancestry, quantity): mean over genes of the predicted correlation sqrt(Var_i(X mu))
	# or of the observed correlation, with its standard error across genes
	rows = list()
	for (target_tissue in names(tissue_colors)) {
		tissue_df = merged_df[merged_df$target_tissue == target_tissue, ]
		if (nrow(tissue_df) < 2) {
			print(paste("WARNING: fewer than 2 genes compared in", target_tissue, "; skipping"))
			next
		}
		for (population in population_levels) {
			suffix = paste0("_", tolower(population))
			for (quantity in names(quantity_colors)) {
				if (quantity == "Predicted r") {
					vals = tissue_df[[paste0("predicted_r", suffix)]]
				} else {
					vals = tissue_df[[paste0("observed_r", suffix)]]
				}
				rows[[length(rows) + 1]] = data.frame(
					target_tissue=target_tissue,
					population=population,
					quantity=quantity,
					n_genes=length(vals),
					mean_r=mean(vals),
					mean_r_se=sd(vals)/sqrt(length(vals)),
					stringsAsFactors=FALSE
				)
			}
		}
	}
	if (length(rows) == 0) {
		return(data.frame())
	}
	summary_df = do.call(rbind, rows)
	summary_df$population = factor(summary_df$population, levels=population_levels)
	summary_df$quantity = factor(summary_df$quantity, levels=names(quantity_colors))
	present_tissues = names(tissue_colors)[names(tissue_colors) %in% summary_df$target_tissue]
	summary_df$tissue = factor(summary_df$target_tissue, levels=present_tissues, labels=tissue_labels[present_tissues])
	return(summary_df)
}


make_population_mean_r_panel_plot <- function(summary_df, ylab) {
	# One panel per tissue: mean predicted and mean observed correlation (points with 95% CI across genes)
	# in Europeans and in African-ancestry individuals, joined by a line per quantity. If the predictions
	# port as expected, the observed line drops (or rises) in parallel with the predicted line.
	quantity_offsets = c("Predicted r"=-0.08, "Observed r"=0.08)
	summary_df$x = as.integer(summary_df$population) + quantity_offsets[as.character(summary_df$quantity)]
	n_gene_df = aggregate(n_genes ~ tissue, data=summary_df, FUN=max)
	n_gene_df$n_label = paste0("N=", n_gene_df$n_genes, " genes")

	pp = ggplot(summary_df, aes(x=x, y=mean_r, color=quantity, group=quantity)) +
		geom_hline(yintercept=0, color="grey80", linewidth=0.4) +
		geom_line(linewidth=0.6) +
		geom_errorbar(aes(ymin=mean_r - 1.96*mean_r_se, ymax=mean_r + 1.96*mean_r_se), width=0.07, linewidth=0.5) +
		geom_point(size=2.0) +
		geom_text(data=n_gene_df, aes(x=Inf, y=Inf, label=n_label), color=reference_color, size=2.6, hjust=1.1, vjust=1.6, inherit.aes=FALSE) +
		facet_wrap(~tissue, nrow=1) +
		scale_x_continuous(breaks=seq_along(population_levels), labels=population_levels, limits=c(0.6, length(population_levels) + 0.4)) +
		scale_color_manual(values=quantity_colors, name=NULL) +
		figure_theme() +
		theme(
			legend.position="bottom",
			legend.key=element_blank(),
			legend.margin=margin(0, 0, 0, 0),
			strip.background=element_blank(),
			strip.text=element_text(face="bold", size=9.5),
			panel.spacing=unit(0.8, "lines")
		) +
		labs(x="Evaluation ancestry", y=ylab)
	return(pp)
}


make_gene_level_eur_vs_afr_scatter_panel_plot <- function(merged_df, eur_col, afr_col, tissue_colors, tissue_labels, xlab, ylab) {
	# One panel per tissue: per-gene value in African-ancestry individuals (y) against Europeans (x), with the
	# y = x line. Each panel is annotated with the gene count and the Pearson correlation across genes.
	plot_df = merged_df[merged_df$target_tissue %in% names(tissue_colors), ]
	if (nrow(plot_df) == 0) {
		return(NULL)
	}
	present_tissues = names(tissue_colors)[names(tissue_colors) %in% plot_df$target_tissue]
	plot_df$tissue = factor(plot_df$target_tissue, levels=present_tissues, labels=tissue_labels[present_tissues])
	plot_df$x_val = plot_df[[eur_col]]
	plot_df$y_val = plot_df[[afr_col]]
	tissue_colors_use = setNames(as.character(tissue_colors[present_tissues]), tissue_labels[present_tissues])

	label_rows = list()
	for (tissue_label in levels(plot_df$tissue)) {
		tissue_df = plot_df[plot_df$tissue == tissue_label, ]
		pearson_r = ifelse(nrow(tissue_df) >= 3, cor(tissue_df$x_val, tissue_df$y_val), NA)
		label_rows[[length(label_rows) + 1]] = data.frame(tissue=tissue_label, label=paste0("N=", nrow(tissue_df), "; r=", sprintf("%.2f", pearson_r)), stringsAsFactors=FALSE)
	}
	label_df = do.call(rbind, label_rows)
	label_df$tissue = factor(label_df$tissue, levels=levels(plot_df$tissue))

	pp = ggplot(plot_df, aes(x=x_val, y=y_val, color=tissue)) +
		geom_point(alpha=0.3, size=0.7, stroke=0) +
		geom_abline(slope=1, intercept=0, linetype="dashed", color=reference_color, linewidth=0.5) +
		geom_text(data=label_df, aes(x=-Inf, y=Inf, label=label), color=reference_color, size=2.6, hjust=-0.1, vjust=1.6, inherit.aes=FALSE) +
		facet_wrap(~tissue, nrow=1) +
		scale_color_manual(values=tissue_colors_use, guide="none") +
		figure_theme() +
		theme(
			strip.background=element_blank(),
			strip.text=element_text(face="bold", size=9.5),
			panel.spacing=unit(0.8, "lines"),
			axis.text=element_text(size=9)
		) +
		labs(x=xlab, y=ylab)
	return(pp)
}


compute_portability_df <- function(merged_df, tissue_colors) {
	# Per tissue: AFR / EUR ratio of the mean predicted correlation and of the mean observed correlation,
	# with delta-method standard errors that account for the pairing of genes across ancestries:
	# Var(a/e) ~ (a/e)^2 (Var(a)/a^2 + Var(e)/e^2 - 2 Cov(a, e)/(a e)), with Var and Cov of the means.
	rows = list()
	for (target_tissue in names(tissue_colors)) {
		tissue_df = merged_df[merged_df$target_tissue == target_tissue, ]
		n_genes = nrow(tissue_df)
		if (n_genes < 3) {
			next
		}
		ratio_values = c()
		for (quantity in c("predicted", "observed")) {
			afr_vals = tissue_df[[paste0(quantity, "_r_afr")]]
			eur_vals = tissue_df[[paste0(quantity, "_r_eur")]]
			afr_mean = mean(afr_vals)
			eur_mean = mean(eur_vals)
			ratio = afr_mean/eur_mean
			rel_var = var(afr_vals)/n_genes/afr_mean^2 + var(eur_vals)/n_genes/eur_mean^2 - 2.0*cov(afr_vals, eur_vals)/n_genes/(afr_mean*eur_mean)
			ratio_values[paste0(quantity, "_ratio")] = ratio
			ratio_values[paste0(quantity, "_ratio_se")] = abs(ratio)*sqrt(max(rel_var, 0.0))
		}
		rows[[length(rows) + 1]] = data.frame(
			target_tissue=target_tissue,
			n_genes=n_genes,
			predicted_ratio=ratio_values[["predicted_ratio"]],
			predicted_ratio_se=ratio_values[["predicted_ratio_se"]],
			observed_ratio=ratio_values[["observed_ratio"]],
			observed_ratio_se=ratio_values[["observed_ratio_se"]],
			stringsAsFactors=FALSE
		)
	}
	if (length(rows) == 0) {
		return(data.frame())
	}
	portability_df = do.call(rbind, rows)
	present_tissues = names(tissue_colors)[names(tissue_colors) %in% portability_df$target_tissue]
	portability_df$tissue = factor(portability_df$target_tissue, levels=present_tissues, labels=gsub("_", " ", present_tissues))
	return(portability_df)
}


make_portability_plot <- function(portability_df, tissue_colors) {
	# One point per tissue: observed AFR/EUR ratio of the mean correlation (y) against the predicted ratio
	# (x), with 95% CIs in both directions and the y = x line (predictions port as expected)
	tissue_colors_use = setNames(as.character(tissue_colors), gsub("_", " ", names(tissue_colors)))
	# Shared x / y limits covering every CI and the point (1, 1)
	axis_values = c(portability_df$observed_ratio - 1.96*portability_df$observed_ratio_se, portability_df$observed_ratio + 1.96*portability_df$observed_ratio_se, portability_df$predicted_ratio - 1.96*portability_df$predicted_ratio_se, portability_df$predicted_ratio + 1.96*portability_df$predicted_ratio_se, 1.0)
	axis_limits = range(axis_values[is.finite(axis_values)])

	pp = ggplot(portability_df, aes(x=predicted_ratio, y=observed_ratio, color=tissue)) +
		geom_abline(slope=1, intercept=0, linetype="dashed", color=reference_color, linewidth=0.5) +
		geom_hline(yintercept=1, color="grey85", linewidth=0.4) +
		geom_vline(xintercept=1, color="grey85", linewidth=0.4) +
		geom_segment(aes(x=predicted_ratio, xend=predicted_ratio, y=observed_ratio - 1.96*observed_ratio_se, yend=observed_ratio + 1.96*observed_ratio_se), linewidth=0.5) +
		geom_segment(aes(x=predicted_ratio - 1.96*predicted_ratio_se, xend=predicted_ratio + 1.96*predicted_ratio_se, y=observed_ratio, yend=observed_ratio), linewidth=0.5) +
		geom_point(size=2.4) +
		scale_color_manual(values=tissue_colors_use, name="Tissue") +
		coord_equal(xlim=axis_limits, ylim=axis_limits) +
		figure_theme() +
		theme(legend.key=element_blank()) +
		labs(x="Predicted AFR / EUR mean r", y="Observed AFR / EUR mean r")
	return(pp)
}





compute_top_n_eur_vs_afr_df <- function(merged_df, model_suffix, min_genes, tissue_colors, tissue_labels) {
	# Top-N-genes curves for the "Predicted r with genetic expression" ranking, once per ancestry: within each
	# tissue, genes (compared in both ancestries) are sorted by decreasing predicted corr(X mu, X beta)
	# = sqrt(Var_i(X mu) / predicted cis-SNP h2), computed in that ancestry's own samples, and the running mean
	# of that ancestry's observed expression correlation over the top N genes is recorded for N >= min_genes.
	# Same score as "Predicted r with genetic expression" in get_gene_ranking_scores of
	# visualize_personalized_expression_correlations.R. A gene missing the score in one ancestry is dropped
	# from that ancestry's curve only. One row per (tissue, ancestry, N).
	tissue_arr = c()
	population_arr = c()
	top_n_arr = c()
	mean_r_arr = c()
	for (target_tissue in names(tissue_colors)) {
		tissue_df = merged_df[merged_df$target_tissue == target_tissue, ]
		for (population in population_levels) {
			suffix = paste0("_", tolower(population))
			score = sqrt(tissue_df[[paste0("rescaled_predicted_expression_variance", suffix)]]/tissue_df[[paste0("predicted_cis_snp_h2", model_suffix, suffix)]])
			observed_r = tissue_df[[paste0("observed_r", suffix)]]
			valid = is.finite(score) & is.finite(observed_r)
			n_valid = sum(valid)
			if (n_valid < min_genes) {
				print(paste("WARNING:", target_tissue, population, ": fewer than", min_genes, "genes with a predicted r with genetic expression; skipping"))
				next
			}
			ordering = order(score[valid], decreasing=TRUE)
			running_mean = cumsum(observed_r[valid][ordering])/seq_len(n_valid)
			top_n = seq_len(n_valid)
			keep = top_n >= min_genes
			tissue_arr = c(tissue_arr, rep(target_tissue, sum(keep)))
			population_arr = c(population_arr, rep(population, sum(keep)))
			top_n_arr = c(top_n_arr, top_n[keep])
			mean_r_arr = c(mean_r_arr, running_mean[keep])
		}
	}
	if (length(tissue_arr) == 0) {
		return(data.frame())
	}
	summary_df = data.frame(
		target_tissue=as.character(tissue_arr),
		population=factor(population_arr, levels=population_levels),
		top_n=as.integer(top_n_arr),
		mean_observed_r=as.numeric(mean_r_arr),
		stringsAsFactors=FALSE
	)
	present_tissues = names(tissue_colors)[names(tissue_colors) %in% summary_df$target_tissue]
	summary_df$tissue = factor(summary_df$target_tissue, levels=present_tissues, labels=tissue_labels[present_tissues])
	return(summary_df)
}


make_top_n_eur_vs_afr_panel_plot <- function(summary_df, xlab, ylab) {
	# Top-N-genes plot styled like make_per_tissue_top_n_mean_observed_r_panel_plot in
	# visualize_personalized_expression_correlations.R: one panel per tissue, mean observed expression
	# correlation of the top N genes (y) against N on a log scale (x), one line per ancestry
	n_gene_df = aggregate(top_n ~ tissue, data=summary_df, FUN=max)
	n_gene_df$n_label = paste0("N=", n_gene_df$top_n, " genes")

	pp = ggplot(summary_df, aes(x=top_n, y=mean_observed_r, color=population)) +
		geom_line(linewidth=0.85, na.rm=TRUE) +
		geom_text(data=n_gene_df, aes(x=Inf, y=Inf, label=n_label), color=reference_color, size=2.6, hjust=1.1, vjust=1.6, inherit.aes=FALSE) +
		facet_wrap(~tissue, nrow=1) +
		scale_color_manual(values=population_colors, labels=population_legend_labels, name="Evaluation ancestry") +
		scale_x_log10(breaks=c(10, 100, 1000, 10000), labels=c("10", "100", "1,000", "10,000")) +
		figure_theme() +
		theme(
			legend.position="bottom",
			legend.key=element_blank(),
			legend.margin=margin(0, 0, 0, 0),
			strip.background=element_blank(),
			strip.text=element_text(face="bold", size=9.5),
			panel.spacing=unit(0.8, "lines"),
			axis.text.x=element_text(size=9)
		) +
		labs(x=xlab, y=ylab)
	return(pp)
}





compute_fsr_bin_summary_df <- function(df, fsr_col, tissue_colors) {
	# One row per (tissue, expression-FSR bin): observed FSR (fraction of genes whose rescaled
	# expression correlation is negative) with a 95% CI, expected FSR (mean of the per-gene FSR in the bin)
	# and the gene count. Same definitions as compute_fsr_bin_summary_df in
	# visualize_personalized_expression_correlations.R (mean-correlation columns omitted).
	df = df[df$target_tissue %in% names(tissue_colors), ]
	df = df[!is.na(df[[fsr_col]]) & !is.na(df$rescaled_expression_correlation), ]
	df$fsr_bin = cut(df[[fsr_col]], breaks=fsr_bin_breaks, labels=fsr_bin_labels, include.lowest=TRUE, right=TRUE)

	rows = list()
	for (target_tissue in names(tissue_colors)) {
		for (bin_iter in seq_along(fsr_bin_labels)) {
			bin_name = fsr_bin_labels[bin_iter]
			bin_df = df[df$target_tissue == target_tissue & df$fsr_bin == bin_name, ]
			n_genes = nrow(bin_df)
			observed_fsr = NA
			observed_fsr_lb = NA
			observed_fsr_ub = NA
			expected_fsr = NA
			if (n_genes > 0) {
				observed_fsr = sum(bin_df$rescaled_expression_correlation < 0)/n_genes
				observed_fsr_se = sqrt(observed_fsr*(1.0 - observed_fsr)/n_genes)
				observed_fsr_lb = max(0.0, observed_fsr - 1.96*observed_fsr_se)
				observed_fsr_ub = min(1.0, observed_fsr + 1.96*observed_fsr_se)
				expected_fsr = mean(bin_df[[fsr_col]])
			}
			rows[[length(rows) + 1]] = data.frame(
				target_tissue=target_tissue,
				fsr_bin=bin_name,
				x=bin_iter,
				n_genes=n_genes,
				observed_fsr=observed_fsr,
				observed_fsr_lb=observed_fsr_lb,
				observed_fsr_ub=observed_fsr_ub,
				expected_fsr=expected_fsr,
				stringsAsFactors=FALSE
			)
		}
	}
	summary_df = do.call(rbind, rows)
	summary_df$fsr_bin = factor(summary_df$fsr_bin, levels=fsr_bin_labels)
	return(summary_df)
}


make_two_row_fsr_calibration_panel_plot <- function(summary_df, tissue_colors, tissue_labels, xlab) {
	# FSR calibration plot styled like make_per_tissue_fsr_calibration_panel_plot in
	# visualize_personalized_expression_correlations.R, with one row per ancestry (European on top, African
	# below) and one column per tissue: observed FSR (points with 95% CI, joined by a line, colored by
	# tissue) against the expected FSR (dashed segment) in each expression-FSR bin, with per-bin gene counts.
	summary_df$tissue = factor(summary_df$target_tissue, levels=names(tissue_colors), labels=tissue_labels[names(tissue_colors)])
	summary_df$population = factor(summary_df$population, levels=population_levels, labels=population_legend_labels[population_levels])
	tissue_colors_use = setNames(as.character(tissue_colors), tissue_labels[names(tissue_colors)])

	# (Per-bin gene counts are not drawn; they are in the bin summary file)
	pp = ggplot(summary_df, aes(x=x, y=observed_fsr, color=tissue)) +
		geom_segment(aes(x=x - 0.35, xend=x + 0.35, y=expected_fsr, yend=expected_fsr), color=reference_color, linetype="dashed", linewidth=0.6, na.rm=TRUE) +
		geom_line(linewidth=0.7, na.rm=TRUE) +
		geom_pointrange(aes(ymin=observed_fsr_lb, ymax=observed_fsr_ub), linewidth=0.7, na.rm=TRUE) +
		geom_point(size=2.4, na.rm=TRUE) +
		facet_grid(population ~ tissue) +
		scale_color_manual(values=tissue_colors_use) +
		scale_x_continuous(breaks=seq_along(fsr_bin_labels), labels=fsr_bin_labels) +
		scale_y_continuous(breaks=seq(0, 0.6, by=0.1)) +
		coord_cartesian(ylim=c(0, 0.6)) +
		figure_theme() +
		theme(
			legend.position="none",
			strip.background=element_blank(),
			strip.text.x=element_text(face="bold", size=9.5),
			strip.text.y=element_text(face="bold", size=11),
			axis.text.x=element_text(angle=35, hjust=1, vjust=1, size=8),
			panel.spacing=unit(0.8, "lines")
		) +
		labs(x=xlab, y="Observed FSR")
	return(pp)
}





compute_predicted_genetic_r_bin_summary_df <- function(df, predicted_r_col, observed_r_col, observed_h2_col, tissue_colors) {
	# One row per (tissue, predefined predicted-correlation bin). Same definitions as
	# compute_predicted_r_bin_summary_df (predefined bins) + compute_predicted_r_bin_stats in
	# visualize_personalized_expression_correlations.R: the bin's mean predicted correlation, and its observed
	# correlation mean(observed r) / sqrt(mean(observed h2)) with a delta-method SE including the
	# numerator-denominator covariance; a bin whose mean observed h2 is not positive gets NA. Every bin is
	# kept (an empty one with NA statistics).
	df = df[df$target_tissue %in% names(tissue_colors), ]
	df = df[!is.na(df[[predicted_r_col]]) & !is.na(df[[observed_r_col]]) & !is.na(df[[observed_h2_col]]), ]

	rows = list()
	for (target_tissue in names(tissue_colors)) {
		tissue_df = df[df$target_tissue == target_tissue, ]
		if (nrow(tissue_df) == 0) {
			next
		}
		tissue_df$r_bin = as.integer(cut(tissue_df[[predicted_r_col]], breaks=predicted_r_bin_breaks, labels=FALSE, include.lowest=TRUE, right=TRUE))
		for (bin_iter in seq_along(predicted_r_bin_labels)) {
			bin_df = tissue_df[!is.na(tissue_df$r_bin) & tissue_df$r_bin == bin_iter, ]
			n_genes = nrow(bin_df)
			mean_predicted = NA
			mean_predicted_se = NA
			mean_observed = NA
			mean_observed_se = NA
			if (n_genes > 0) {
				predicted_r = bin_df[[predicted_r_col]]
				observed_r = bin_df[[observed_r_col]]
				observed_h2 = bin_df[[observed_h2_col]]
				mean_predicted = mean(predicted_r)
				mean_r = mean(observed_r)
				mean_h2 = mean(observed_h2)
				if (mean_h2 > 0) {
					mean_observed = mean_r/sqrt(mean_h2)
				} else {
					print(paste0("WARNING: ", target_tissue, " ", predicted_r_col, " bin ", bin_iter, ": mean ", observed_h2_col, " = ", mean_h2, " is not positive; observed r / sqrt(h2) set to NA"))
				}
				if (n_genes > 1) {
					mean_predicted_se = sd(predicted_r)/sqrt(n_genes)
					if (mean_h2 > 0) {
						var_mean_r = var(observed_r)/n_genes
						var_mean_h2 = var(observed_h2)/n_genes
						cov_mean_r_h2 = cov(observed_r, observed_h2)/n_genes
						# Delta method for g(a, b) = a / sqrt(b): dg/da = 1/sqrt(b), dg/db = -a / (2 b^(3/2))
						var_ratio = var_mean_r/mean_h2 - mean_r*cov_mean_r_h2/mean_h2^2 + mean_r^2*var_mean_h2/(4.0*mean_h2^3)
						mean_observed_se = sqrt(max(var_ratio, 0.0))
					}
				}
			}
			rows[[length(rows) + 1]] = data.frame(
				target_tissue=target_tissue,
				predicted_r_bin=bin_iter,
				predicted_r_bin_label=predicted_r_bin_labels[bin_iter],
				n_genes=n_genes,
				mean_predicted_r=mean_predicted,
				mean_predicted_r_lb=mean_predicted - 1.96*mean_predicted_se,
				mean_predicted_r_ub=mean_predicted + 1.96*mean_predicted_se,
				mean_observed_r=mean_observed,
				mean_observed_r_lb=mean_observed - 1.96*mean_observed_se,
				mean_observed_r_ub=mean_observed + 1.96*mean_observed_se,
				stringsAsFactors=FALSE
			)
		}
	}
	if (length(rows) == 0) {
		return(data.frame())
	}
	return(do.call(rbind, rows))
}


make_two_row_predicted_r_binned_panel_plot <- function(summary_df, tissue_colors, tissue_labels, xlab, ylab, y_max_cap) {
	# Predicted r calibration plot styled like make_per_tissue_predicted_vs_observed_r_binned_panel_plot in
	# visualize_personalized_expression_correlations.R, with one row per ancestry (European on top, African
	# below) and one column per tissue: observed r / sqrt(h2) (points with 95% CI, joined by a line, colored
	# by tissue) against the bin's mean predicted correlation (dashed segment) in each predefined bin.
	# The y-axis is shared by all panels and fit to the data, capped at y_max_cap.
	summary_df$tissue = factor(summary_df$target_tissue, levels=names(tissue_colors), labels=tissue_labels[names(tissue_colors)])
	summary_df$population = factor(summary_df$population, levels=population_levels, labels=population_legend_labels[population_levels])
	tissue_colors_use = setNames(as.character(tissue_colors), tissue_labels[names(tissue_colors)])

	y_max = ceiling(max(c(summary_df$mean_observed_r_ub, summary_df$mean_observed_r, summary_df$mean_predicted_r), na.rm=TRUE)/0.05)*0.05
	if (y_max > y_max_cap) {
		print(paste("Two-row predicted r calibration plot: y-axis capped at", y_max_cap, "(data reach", y_max, ")"))
		y_max = y_max_cap
	}
	y_breaks = pretty(c(0, y_max), n=4)

	pp = ggplot(summary_df, aes(x=predicted_r_bin, y=mean_observed_r, color=tissue)) +
		geom_segment(aes(x=predicted_r_bin - 0.35, xend=predicted_r_bin + 0.35, y=mean_predicted_r, yend=mean_predicted_r), color=reference_color, linetype="dashed", linewidth=0.6, na.rm=TRUE) +
		geom_line(linewidth=0.7, na.rm=TRUE) +
		geom_pointrange(aes(ymin=mean_observed_r_lb, ymax=mean_observed_r_ub), linewidth=0.7, na.rm=TRUE) +
		geom_point(size=2.4, na.rm=TRUE) +
		facet_grid(population ~ tissue) +
		scale_color_manual(values=tissue_colors_use) +
		scale_x_continuous(breaks=seq_along(predicted_r_bin_labels), labels=predicted_r_bin_labels) +
		scale_y_continuous(breaks=y_breaks, expand=expansion(mult=c(0.02, 0.03))) +
		coord_cartesian(ylim=c(0, y_max)) +
		figure_theme() +
		theme(
			legend.position="none",
			strip.background=element_blank(),
			strip.text.x=element_text(face="bold", size=9.5),
			strip.text.y=element_text(face="bold", size=11),
			axis.text.x=element_text(angle=35, hjust=1, vjust=1, size=8),
			panel.spacing=unit(0.8, "lines")
		) +
		labs(x=xlab, y=ylab)
	return(pp)
}





compute_paired_eur_vs_afr_test_df <- function(merged_df, eur_col, afr_col, tissue_colors) {
	# Per tissue, over genes with finite values in both ancestries: gene count, mean and median in each
	# ancestry, and the two-sided paired Wilcoxon signed-rank test of EUR vs AFR (same genes in both)
	rows = list()
	for (target_tissue in names(tissue_colors)) {
		tissue_df = merged_df[merged_df$target_tissue == target_tissue, ]
		valid = is.finite(tissue_df[[eur_col]]) & is.finite(tissue_df[[afr_col]])
		eur_vals = tissue_df[[eur_col]][valid]
		afr_vals = tissue_df[[afr_col]][valid]
		pvalue = NA
		if (sum(valid) >= 3) {
			pvalue = wilcox.test(eur_vals, afr_vals, paired=TRUE, exact=FALSE)$p.value
		}
		rows[[length(rows) + 1]] = data.frame(
			target_tissue=target_tissue,
			n_genes=sum(valid),
			mean_eur=ifelse(sum(valid) > 0, mean(eur_vals), NA),
			mean_afr=ifelse(sum(valid) > 0, mean(afr_vals), NA),
			median_eur=ifelse(sum(valid) > 0, median(eur_vals), NA),
			median_afr=ifelse(sum(valid) > 0, median(afr_vals), NA),
			paired_wilcoxon_pvalue=pvalue,
			stringsAsFactors=FALSE
		)
	}
	return(do.call(rbind, rows))
}


make_two_row_eur_vs_afr_histogram_plot <- function(merged_df, eur_col, afr_col, test_df, tissue_colors, tissue_labels, xlab, binwidth) {
	# Histograms of a per-gene quantity with one row per ancestry (European on top, African below) and one
	# column per tissue, over genes with finite values in both ancestries (same genes in both rows). The dashed
	# line marks each panel's median. Column titles carry the tissue, sample sizes and the paired Wilcoxon p-value.
	long_rows = list()
	for (population in population_levels) {
		value_col = ifelse(population == "EUR", eur_col, afr_col)
		valid = is.finite(merged_df[[eur_col]]) & is.finite(merged_df[[afr_col]]) & merged_df$target_tissue %in% names(tissue_colors)
		long_rows[[population]] = data.frame(target_tissue=merged_df$target_tissue[valid], population=population, value=merged_df[[value_col]][valid], stringsAsFactors=FALSE)
	}
	plot_df = do.call(rbind, long_rows)
	if (nrow(plot_df) == 0) {
		return(NULL)
	}

	# Column titles: tissue (with sample sizes) and the paired test p-value
	column_labels = c()
	for (target_tissue in names(tissue_colors)) {
		pvalue = test_df$paired_wilcoxon_pvalue[test_df$target_tissue == target_tissue]
		pvalue_str = ifelse(is.na(pvalue), "NA", format.pval(pvalue, digits=2, eps=1e-300))
		column_labels[target_tissue] = paste0(tissue_labels[[target_tissue]], "\npaired Wilcoxon p = ", pvalue_str)
	}
	present_tissues = names(tissue_colors)[names(tissue_colors) %in% plot_df$target_tissue]
	plot_df$tissue = factor(plot_df$target_tissue, levels=present_tissues, labels=column_labels[present_tissues])
	plot_df$population = factor(plot_df$population, levels=population_levels, labels=population_legend_labels[population_levels])
	median_df = aggregate(value ~ tissue + population, data=plot_df, FUN=median)
	tissue_colors_use = setNames(as.character(tissue_colors[present_tissues]), column_labels[present_tissues])

	pp = ggplot(plot_df, aes(x=value, fill=tissue)) +
		geom_histogram(binwidth=binwidth, boundary=0, color=NA) +
		geom_vline(data=median_df, aes(xintercept=value), color=reference_color, linetype="dashed", linewidth=0.5) +
		facet_grid(population ~ tissue) +
		scale_fill_manual(values=tissue_colors_use, guide="none") +
		figure_theme() +
		theme(
			strip.background=element_blank(),
			strip.text.x=element_text(face="bold", size=9),
			strip.text.y=element_text(face="bold", size=11),
			axis.text=element_text(size=9),
			panel.spacing=unit(0.8, "lines")
		) +
		labs(x=xlab, y="Number of genes")
	return(pp)
}





#####################
# Command line args
#####################
eur_per_tissue_expression_correlation_dir = args[1]
afr_per_tissue_expression_correlation_dir = args[2]
tissue_names_file = args[3]
visualization_dir = args[4]


#####################
# Load in data
#####################
tissue_info_df = read.table(tissue_names_file, header=TRUE, sep="\t")

eur_results_df = load_per_tissue_expression_correlations(eur_per_tissue_expression_correlation_dir, tissue_info_df)
afr_results_df = load_per_tissue_expression_correlations(afr_per_tissue_expression_correlation_dir, tissue_info_df)

merged_df = merge_eur_afr_results(eur_results_df, afr_results_df)


#####################
# Five tissues used for plotting (same colors as visualize_personalized_expression_correlations.R)
#####################
five_tissue_colors = c(
	"Heart_Left_Ventricle"="#4C72B0",
	"Brain_Cortex"="#DD8452",
	"Liver"="#55A868",
	"Whole_Blood"="#C44E52",
	"Muscle_Skeletal"="#8172B3"
)
missing_tissues = setdiff(names(five_tissue_colors), unique(merged_df$target_tissue))
if (length(missing_tissues) > 0) {
	print(paste("WARNING: no genes compared in both ancestries for:", paste(missing_tissues, collapse=", ")))
}

sample_size_df = compute_tissue_sample_size_df(eur_results_df, afr_results_df, merged_df, five_tissue_colors)
print(sample_size_df)
write.table(sample_size_df, paste0(visualization_dir, "five_tissue_eur_vs_afr_sample_sizes.txt"), quote=FALSE, sep="\t", row.names=FALSE)
tissue_labels = make_tissue_labels(sample_size_df, five_tissue_colors)


#####################
# Two-row expression-FSR (AF-specific) calibration plot: European row restricted to genes heritable in
# Europeans, African row restricted to genes heritable in African-ancestry individuals (each ancestry's own
# cis-h2 LRT p < threshold, over all genes analyzed in that ancestry, not only genes shared by both)
#####################
two_row_fsr_col = "expression_FSR_af_specific"
two_row_fsr_summary_df = data.frame()
for (population in population_levels) {
	if (population == "EUR") {
		population_df = eur_results_df
	} else {
		population_df = afr_results_df
	}
	population_heritable_df = population_df[!is.na(population_df$cis_snp_h2_pvalue) & population_df$cis_snp_h2_pvalue < heritability_pvalue_threshold, ]
	print(paste(population, "genes heritable in", population, "(p <", heritability_pvalue_threshold, "):", nrow(population_heritable_df), "of", nrow(population_df)))
	population_summary_df = compute_fsr_bin_summary_df(population_heritable_df, two_row_fsr_col, five_tissue_colors)
	population_summary_df$population = population
	two_row_fsr_summary_df = rbind(two_row_fsr_summary_df, population_summary_df)
}
print(two_row_fsr_summary_df)
write.table(two_row_fsr_summary_df, paste0(visualization_dir, "five_tissue_eur_vs_afr_", two_row_fsr_col, "_bin_summary_own_ancestry_heritable.txt"), quote=FALSE, sep="\t", row.names=FALSE)
two_row_fsr_plot = make_two_row_fsr_calibration_panel_plot(two_row_fsr_summary_df, five_tissue_colors, tissue_labels, "Expression-FSR (AF-specific)")
ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_", two_row_fsr_col, "_calibration_per_tissue_panels_own_ancestry_heritable.pdf"), two_row_fsr_plot, width=10.0, height=5.6)


#####################
# Two-row histograms (European top, African bottom) of per-gene expected quantities under the AF-specific
# model, over all genes compared in both ancestries, with a paired Wilcoxon signed-rank test per tissue:
#   (1) expected FSR (expression_FSR_af_specific)
#   (2) expected correlation with genetic expression (predicted_genetic_r_af_specific)
#####################
histogram_definitions = list(
	expression_FSR_af_specific=list(xlab="Expected FSR (AF-specific)", binwidth=0.01),
	predicted_genetic_r_af_specific=list(xlab="Expected r with genetic expression (AF-specific)", binwidth=0.02)
)
for (histogram_col in names(histogram_definitions)) {
	eur_col = paste0(histogram_col, "_eur")
	afr_col = paste0(histogram_col, "_afr")
	print(paste("Plotting EUR vs AFR histogram of", histogram_col, "for all genes compared in both ancestries"))
	test_df = compute_paired_eur_vs_afr_test_df(merged_df, eur_col, afr_col, five_tissue_colors)
	print(test_df)
	write.table(test_df, paste0(visualization_dir, "five_tissue_eur_vs_afr_", histogram_col, "_paired_test_all_genes.txt"), quote=FALSE, sep="\t", row.names=FALSE)
	histogram_plot = make_two_row_eur_vs_afr_histogram_plot(merged_df, eur_col, afr_col, test_df, five_tissue_colors, tissue_labels, histogram_definitions[[histogram_col]]$xlab, histogram_definitions[[histogram_col]]$binwidth)
	if (is.null(histogram_plot)) {
		print(paste("WARNING: no genes with", histogram_col, "in both ancestries; skipping the histogram"))
	} else {
		ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_", histogram_col, "_histogram_all_genes.pdf"), histogram_plot, width=10.5, height=5.6)
	}
}


#####################
# Two-row predicted r with genetic expression (AF-specific) calibration plot, predefined bins, all genes:
# European row from the European results, African row from the African results; within each ancestry the
# predicted r = sqrt(Var_i(X mu) / predicted cis-SNP h2 (AF-specific)) and the observed r / sqrt(observed h2)
# (Haseman-Elston h2) are both computed in that ancestry's own samples. All genes analyzed in each ancestry.
#####################
two_row_predicted_r_col = "predicted_genetic_r_af_specific"
two_row_predicted_r_summary_df = data.frame()
for (population in population_levels) {
	if (population == "EUR") {
		population_df = eur_results_df
	} else {
		population_df = afr_results_df
	}
	population_df[[two_row_predicted_r_col]] = sqrt(population_df$rescaled_predicted_expression_variance/population_df$predicted_cis_snp_h2_af_specific)
	population_summary_df = compute_predicted_genetic_r_bin_summary_df(population_df, two_row_predicted_r_col, observed_r_col, observed_h2_col, five_tissue_colors)
	if (nrow(population_summary_df) == 0) {
		print(paste("WARNING: no", population, "genes for the two-row predicted r calibration plot"))
		next
	}
	population_summary_df$population = population
	two_row_predicted_r_summary_df = rbind(two_row_predicted_r_summary_df, population_summary_df)
}
print(two_row_predicted_r_summary_df)
if (nrow(two_row_predicted_r_summary_df) > 0) {
	write.table(two_row_predicted_r_summary_df, paste0(visualization_dir, "five_tissue_eur_vs_afr_", two_row_predicted_r_col, "_predefined_bin_summary_all_genes.txt"), quote=FALSE, sep="\t", row.names=FALSE)
	two_row_predicted_r_plot = make_two_row_predicted_r_binned_panel_plot(two_row_predicted_r_summary_df, five_tissue_colors, tissue_labels, "Predicted r with genetic expression (AF-specific)", "Observed r / sqrt(h2)", two_row_predicted_r_y_max_cap)
	ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_", two_row_predicted_r_col, "_calibration_per_tissue_panels_predefined_bins_all_genes.pdf"), two_row_predicted_r_plot, width=10.0, height=5.6)
}


#####################
# Plots, once for genes heritable in Europeans (main; no file suffix) and once for all genes compared in
# both ancestries ("_all_genes" suffix)
#####################
heritable_merged_df = merged_df[!is.na(merged_df$cis_snp_h2_pvalue_eur) & merged_df$cis_snp_h2_pvalue_eur < heritability_pvalue_threshold, ]
print(paste("Genes heritable in Europeans (p <", heritability_pvalue_threshold, "):", nrow(heritable_merged_df), "of", nrow(merged_df)))

gene_sets = list(
	heritable_in_eur_genes=list(df=heritable_merged_df, suffix=""),
	all_genes=list(df=merged_df, suffix="_all_genes")
)

for (gene_set_name in names(gene_sets)) {
	gene_set_df = gene_sets[[gene_set_name]]$df
	file_suffix = gene_sets[[gene_set_name]]$suffix
	print(paste("Plotting EUR vs AFR comparisons for", gene_set_name, "(", nrow(gene_set_df), "genes )"))

	# (1) Mean predicted and observed correlation in each ancestry
	mean_r_summary_df = compute_population_mean_r_summary_df(gene_set_df, five_tissue_colors, tissue_labels)
	print(mean_r_summary_df)
	if (nrow(mean_r_summary_df) == 0) {
		print("WARNING: no tissues with genes compared in both ancestries; skipping the mean correlation plot")
	} else {
		write.table(mean_r_summary_df, paste0(visualization_dir, "five_tissue_eur_vs_afr_mean_r_summary", file_suffix, ".txt"), quote=FALSE, sep="\t", row.names=FALSE)
		mean_r_plot = make_population_mean_r_panel_plot(mean_r_summary_df, "Mean r with observed expression")
		ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_mean_r_per_tissue_panels", file_suffix, ".pdf"), mean_r_plot, width=9.5, height=3.3)
	}

	# (2) Per-gene predicted correlation (LD / allele-frequency change only) and observed correlation, AFR vs EUR
	predicted_scatter_plot = make_gene_level_eur_vs_afr_scatter_panel_plot(gene_set_df, "predicted_r_eur", "predicted_r_afr", five_tissue_colors, tissue_labels, "Predicted r (EUR)", "Predicted r (AFR)")
	if (!is.null(predicted_scatter_plot)) {
		ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_predicted_r_gene_scatter", file_suffix, ".pdf"), predicted_scatter_plot, width=10.5, height=2.9)
	}
	observed_scatter_plot = make_gene_level_eur_vs_afr_scatter_panel_plot(gene_set_df, "observed_r_eur", "observed_r_afr", five_tissue_colors, tissue_labels, "Observed r (EUR)", "Observed r (AFR)")
	if (!is.null(observed_scatter_plot)) {
		ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_observed_r_gene_scatter", file_suffix, ".pdf"), observed_scatter_plot, width=10.5, height=2.9)
	}

	# (3) Portability: observed vs predicted AFR / EUR ratio of the mean correlation, one point per tissue
	portability_df = compute_portability_df(gene_set_df, five_tissue_colors)
	print(portability_df)
	if (nrow(portability_df) == 0) {
		print("WARNING: no tissues with at least 3 genes compared in both ancestries; skipping the portability plot")
	} else {
		write.table(portability_df, paste0(visualization_dir, "five_tissue_eur_vs_afr_portability_summary", file_suffix, ".txt"), quote=FALSE, sep="\t", row.names=FALSE)
		portability_plot = make_portability_plot(portability_df, five_tissue_colors)
		ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_portability", file_suffix, ".pdf"), portability_plot, width=5.0, height=4.0)
	}

	# (4) Top-N-genes plot for the "Predicted r with genetic expression" ranking, one line per ancestry,
	# once per residual-variance model
	for (model_suffix in ranking_model_suffixes) {
		print(paste0("Plotting EUR vs AFR top-N-genes mean observed r (model suffix '", model_suffix, "') for ", gene_set_name))
		top_n_df = compute_top_n_eur_vs_afr_df(gene_set_df, model_suffix, top_n_min_genes, five_tissue_colors, tissue_labels)
		if (nrow(top_n_df) == 0) {
			print("WARNING: no genes available for the EUR vs AFR top-N-genes plot; skipping")
			next
		}
		write.table(top_n_df, paste0(visualization_dir, "five_tissue_eur_vs_afr_top_n_genes_mean_observed_r_genetic_r", model_suffix, file_suffix, ".txt"), quote=FALSE, sep="\t", row.names=FALSE)
		top_n_plot = make_top_n_eur_vs_afr_panel_plot(top_n_df, "Top N genes (ranked by predicted r with genetic expression)", "Mean observed r of top N genes")
		ggsave(paste0(visualization_dir, "five_tissue_eur_vs_afr_top_n_genes_mean_observed_r_genetic_r", model_suffix, "_per_tissue_panels", file_suffix, ".pdf"), top_n_plot, width=10.5, height=3.5)
	}
}
