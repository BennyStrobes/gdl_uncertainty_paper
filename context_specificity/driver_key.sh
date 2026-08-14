####################
# Input data
####################

# Directory of expression data
gtex_expr_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/residualized_expression/"

# Directory containing genotype data
processed_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Directory containing eQTL summary statistics
eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/eqtl_results/"

# Directory containing borzoi gtex target indices and names
borzoi_gtex_unique_target_names_file="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/targets_gtex_eqtl_only_unique_ordered.txt"


# SLDMC Output root directory (created from SLDMC analysis code base)
sldmc_root="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/sldmc_analysis/"
borzoi_output_dir=${sldmc_root}"processed_borzoi/"
sldmc_results_output_dir=${sldmc_root}"sldmc_results/"


####################
# Output data
####################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/context_specificity/"

expression_differences_results_dir=${output_root}"expression_differences_results/"

visualize_expression_differences_results_dir=${output_root}"visualize_expression_differences_results/"



####################
# Tissues (with representative borzoi target sample) to compute all pairwise expression differences between
####################
tissue_sample_pairs="Heart_Left_Ventricle:GTEX-18465-0926-SM-731AY.1 Brain_Cortex:GTEX-1H3O1-1726-SM-9WYSR.1 Liver:GTEX-11EQ9-0526-SM-5A5JZ.1 Whole_Blood:GTEX-1LB8K-0005-SM-DIPED.1 Muscle_Skeletal:GTEX-13QJ3-0726-SM-5SI68.1"

anno_method="borzoi_finer_magnitude_bins"

genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
sldmc_results_file=${sldmc_results_output_dir}"sldmc_results_cross_tissue_meta_analyzed_default_bootstrap_stats.txt"


# Loop through all (unordered) pairs of the above tissues
if false; then
ii=0
for tissue_sample1 in $tissue_sample_pairs; do
	ii=$((ii+1))
	jj=0
	for tissue_sample2 in $tissue_sample_pairs; do
		jj=$((jj+1))
		# Only run each pair once
		if [ $jj -le $ii ]; then
			continue
		fi

		target_tissue1=${tissue_sample1%%:*}
		target_sample1=${tissue_sample1#*:}
		target_tissue2=${tissue_sample2%%:*}
		target_sample2=${tissue_sample2#*:}

		borzoi_effect_file1=${borzoi_output_dir}${target_tissue1}"_"${target_sample1}"_borzoi_effects.txt.gz"
		borzoi_effect_file2=${borzoi_output_dir}${target_tissue2}"_"${target_sample2}"_borzoi_effects.txt.gz"
		borzoi_annotation_file1=${borzoi_output_dir}${target_tissue1}"_"${target_sample1}"_annotations_default.txt.gz"
		borzoi_annotation_file2=${borzoi_output_dir}${target_tissue2}"_"${target_sample2}"_annotations_default.txt.gz"

		genotype_sample_mapping_file1=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue1}"_expression_samples.txt"
		genotype_sample_mapping_file2=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue2}"_expression_samples.txt"

		expr_file1=${gtex_expr_dir}${target_tissue1}".v10.residualized_expression_renormalized.bed"
		expr_file2=${gtex_expr_dir}${target_tissue2}".v10.residualized_expression_renormalized.bed"

		expr_differences_output_file=${expression_differences_results_dir}"expr_differences_summary_"${target_tissue1}"_"${target_sample1}"_vs_"${target_tissue2}"_"${target_sample2}"_"${anno_method}".txt.gz"

		sbatch predict_expression_differences_across_tissue_pairs_in_shared_individuals.sh $genotype_stem $sldmc_results_file $expr_differences_output_file $anno_method $borzoi_effect_file1 $borzoi_annotation_file1 $genotype_sample_mapping_file1 $expr_file1 $borzoi_effect_file2 $borzoi_annotation_file2 $genotype_sample_mapping_file2 $expr_file2
	done
done
fi

####################
# Merge all pairwise expression differences results into a summary file
####################
sh merge_expression_differences_results.sh $expression_differences_results_dir $anno_method



####################
# Visualize expression differences results
####################
source ~/.bashrc
conda activate plink_env
Rscript visualize_expression_differences_results.R $expression_differences_results_dir $anno_method $visualize_expression_differences_results_dir

