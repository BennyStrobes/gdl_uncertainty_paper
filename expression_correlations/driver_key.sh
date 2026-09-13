#########################
# Input data
#########################

sldmc_summary_file="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/sldmc_analysis/sldmc_results/sldmc_results_cross_tissue_meta_analyzed_default_bootstrap_stats.txt"

# Following directory contains expression files per tissue like: Whole_Blood.v8.residualized_expression_renormalized.bed
gtex_expression_directory="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/gtex_eqtl_expression_processing/residualized_expression/"

# Following contains genotype data
# Chrom specific plink files: gtex_v9_eqtl_chr10.bed
# And then file containing sample indices: genotype_sample_mapping_to_Whole_Blood_expression_samples.txt
gtex_genotype_directory="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Following directory contains Borzoi predicted effect sizes
# Contains one file for each tissue like: "Whole_Blood_GTEX-1LB8K-0005-SM-DIPED.1_borzoi_effects.txt.gz"
borzoi_predicted_effect_sizes_directory="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/sldmc_analysis/processed_borzoi/"

# Directory containing results of borzoi runs
borzoi_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/"

# Directory containing borzoi gtex target indices and names
borzoi_gtex_unique_target_names_file=${borzoi_results_dir}"targets_gtex_v8_eqtl_only_unique_ordered.txt"


#########################
# Output data
#########################
expression_correlations_root_directory="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/expression_correlations/"

per_tissue_personalized_expression_dir=${expression_correlations_root_directory}"personalized_expression_per_tissue/"


#########################
# Code
#########################


if false; then
tail -n +2 "$borzoi_gtex_unique_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
    borzoi_results_file=${borzoi_predicted_effect_sizes_directory}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
    expression_file=${gtex_expression_directory}${target_tissue}".v8.residualized_expression_renormalized.bed"
    plink_genotype_stem=${gtex_genotype_directory}"gtex_v9_eqtl_chr"
    genotype_sample_mapping_file=${gtex_genotype_directory}"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

    expression_correlation_output_file=${per_tissue_personalized_expression_dir}${target_tissue}"_"${target_sample}"_personalized_expression_prediction.txt"

    sbatch personalized_expression_correlations_per_tissue.sh $sldmc_summary_file $borzoi_results_file $expression_file $plink_genotype_stem $genotype_sample_mapping_file $expression_correlation_output_file
done
fi