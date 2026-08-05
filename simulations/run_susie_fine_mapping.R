suppressPackageStartupMessages(library(susieR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
	stop("Usage: Rscript run_susie_fine_mapping.R X_file y_file variant_file gene_name output_file L")
}

x_file <- args[[1]]
y_file <- args[[2]]
variant_file <- args[[3]]
gene_name <- args[[4]]
output_file <- args[[5]]
num_effects <- as.integer(args[[6]])

X <- as.matrix(read.table(x_file, header = FALSE, sep = "\t"))
y <- scan(y_file, quiet = TRUE)
variant_names <- scan(variant_file, what = "character", quiet = TRUE)

if (nrow(X) != length(y)) {
	stop("X rows must equal length(y)")
}
if (ncol(X) != length(variant_names)) {
	stop("X columns must equal length(variant_names)")
}

fit <- susie(X, y, L = num_effects, standardize = FALSE, intercept = TRUE)

pip <- susie_get_pip(fit)
posterior_mean <- colSums(fit$alpha * fit$mu)
posterior_variance <- colSums((fit$alpha * fit$mu2) - (fit$alpha * fit$mu)^2)
posterior_sd <- sqrt(pmax(posterior_variance, 0))

credible_sets <- rep("NA", length(variant_names))
cs_min_abs_corr <- rep(NA_real_, length(variant_names))

if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
	for (cs_iter in seq_along(fit$sets$cs)) {
		cs_name <- names(fit$sets$cs)[[cs_iter]]
		if (is.null(cs_name) || is.na(cs_name) || cs_name == "") {
			cs_name <- paste0("CS", cs_iter)
		}

		cs_indices <- fit$sets$cs[[cs_iter]]
		old_assignments <- credible_sets[cs_indices]
		credible_sets[cs_indices] <- ifelse(old_assignments == "NA", cs_name, paste(old_assignments, cs_name, sep = ","))

		if (!is.null(fit$sets$purity) && cs_name %in% rownames(fit$sets$purity) && "min.abs.corr" %in% colnames(fit$sets$purity)) {
			cs_min_abs_corr[cs_indices] <- fit$sets$purity[cs_name, "min.abs.corr"]
		}
	}
}

results <- data.frame(
	gene = gene_name,
	variant = variant_names,
	susie_pip = pip,
	susie_posterior_mean = posterior_mean,
	susie_posterior_sd = posterior_sd,
	susie_credible_sets = credible_sets,
	susie_cs_min_abs_corr = cs_min_abs_corr
)

write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
