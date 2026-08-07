#################
# Input data
#################


# Directory containing results of borzoi runs
borzoi_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/"

# Directory containing borzoi gtex target indices and names
borzoi_gtex_unique_target_names_file=${borzoi_results_dir}"targets_gtex_eqtl_only_unique_ordered.txt"

# Directory containing genotype data
processed_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/plink_processed_genotype/"

# Directory of expression data
gtex_expr_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/residualized_expression/"

# Directory containing eQTL summary statistics
eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/gtex_eqtl_expression_processing/eqtl_results/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"

# Simulation results directory (for paper visualization purposes)
simulation_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/simulations/corr_sim_inference_results/"
simulation_oracle_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/simulations/corr_sim_borzoi_eqtl_effects/"

# Baseline LD annotations directory (used for variant annotations)
baselineLD_anno_dir="/lab-share/CHIP-Strober-e2/Public/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"

# Requires list of annotations to run S-LDMC on
annotation_name_file="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/sldmc_analysis/input_data/s_ldmc_annotations.txt"

# Gene set anno file
gene_set_anno_file="/lab-share/CHIP-Strober-e2/Public/gene_set_annotations/non_disease_specific_geneset.csv"

# Directory containing SLDMC code
sldmc_code_dir="/lab-share/CHIP-Strober-e2/Public/ben/SLDMC/"

#################
# Output directories
#################
# Output root directory
output_root="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/sldmc_analysis/"

borzoi_output_dir=${output_root}"processed_borzoi/"

bootstrapped_cross_tissue_gene_sets_dir=${output_root}"bootstrapped_gene_sets/"

sldmc_results_output_dir=${output_root}"sldmc_results/"

tissue_permuted_sldmc_results_output_dir=${output_root}"tissue_permuted_sldmc_results/"

visualize_sldmc_results_dir=${output_root}"visualize_sldmc/"





#################
# Code
#################

#################
# Preprocess data to get into correct format for LD corr
#################

#####
# 1. eQTLs have already been processed

#####
# 2. Borzoi effects
if false; then
tail -n +2 "$borzoi_gtex_unique_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	sbatch preprocess_borzoi_data_for_sldmc.sh $gtex_v10_pc_genes_gtf $borzoi_results_dir $borzoi_target_index $gtex_tissue $target_identifier $borzoi_output_dir
done
fi





#####
# 3. Annotation effects
if false; then
tail -n +2 "$borzoi_gtex_unique_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_annotation_filestem=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_annotations"
	sbatch annotate_variant_gene_pairs.sh $borzoi_effect_file $annotation_name_file $borzoi_annotation_filestem $eqtl_sumstats_file $baselineLD_anno_dir $gene_set_anno_file
done
fi



#################
# 4. Generate cross tissue gene sets (and bootstrapped gene sets)
#################
if false; then
sbatch generate_cross_tissue_bootstrapped_gene_sets.sh ${eqtl_sumstats_dir} ${borzoi_output_dir} ${borzoi_gtex_unique_target_names_file} ${bootstrapped_cross_tissue_gene_sets_dir}
fi


#################
# 5. Run LD-corr
#################
annotation_versions="default magnitude_stratified"
if false; then
tail -n +2 "$borzoi_gtex_unique_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	borzoi_effect_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_borzoi_effects.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"

	for annotation_version in $annotation_versions; do
		borzoi_annotation_file=${borzoi_output_dir}${target_tissue}"_"${target_sample}"_annotations_"${annotation_version}".txt.gz"
		sldmc_output_stem=${sldmc_results_output_dir}"sldmc_results_"${target_tissue}"_"${target_sample}"_"${annotation_version}
		sbatch run_sldmc.sh $borzoi_effect_file $eqtl_sumstats_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file ${bootstrapped_cross_tissue_gene_sets_dir}"cross_tissue_gene_set_bootstrap_" $sldmc_output_stem $sldmc_code_dir
	done
done
fi



#################
# 6. Meta-analyze ld-corr results across tissues
#################
# a. create file with list of all output files (one list per annotation version)
if false; then
for annotation_version in $annotation_versions; do
	sldmc_output_file_list=${sldmc_results_output_dir}"sldmc_per_tissue_output_file_list_"${annotation_version}".txt"
	> "$sldmc_output_file_list"
	tail -n +2 "$borzoi_gtex_unique_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_sample target_description target_tissue; do
		sldmc_output_file=${sldmc_results_output_dir}"sldmc_results_"${target_tissue}"_"${target_sample}"_"${annotation_version}"_bootstrap_summary.txt"
		echo "$sldmc_output_file" >> "$sldmc_output_file_list"
	done
done
fi



# b. Meta-analyze LD-corr results across tissues
if false; then
for annotation_version in $annotation_versions; do
	sldmc_output_file_list=${sldmc_results_output_dir}"sldmc_per_tissue_output_file_list_"${annotation_version}".txt"
	meta_analyzed_sldmc_output_stem=${sldmc_results_output_dir}"sldmc_results_cross_tissue_meta_analyzed_"${annotation_version}
	sh meta_analyze_sldmc_results.sh $sldmc_output_file_list $meta_analyzed_sldmc_output_stem
done
fi


# c. Difference each cell from the intercept (per-tissue and cross-tissue meta-analyzed).
# Run separately per annotation version: the default version differences against a single global
# intercept, the magnitude_stratified version against the per-magnitude-bin intercept.
if false; then
for annotation_version in $annotation_versions; do
	sldmc_output_file_list=${sldmc_results_output_dir}"sldmc_per_tissue_output_file_list_"${annotation_version}".txt"
	intercept_diff_output_stem=${sldmc_results_output_dir}"sldmc_results_cross_tissue_meta_analyzed_"${annotation_version}
	sh compute_intercept_differences.sh $sldmc_output_file_list $intercept_diff_output_stem $annotation_version
done
fi



#################
# 5.Tissue permuted run of S-LDMC
#################


# a. Draw the permuted pairing: one row per tissue, holding that tissue plus a randomly drawn partner
# tissue/sample. Written once and reused, so every downstream permuted run uses the same pairing.
tissue_permuted_pairs_file=${tissue_permuted_sldmc_results_output_dir}"tissue_permuted_pairs.txt"
if false; then
python generate_tissue_permuted_pairs.py \
	--borzoi-gtex-unique-target-names-file $borzoi_gtex_unique_target_names_file \
	--tissue-permuted-pairs-output-file $tissue_permuted_pairs_file \
	--seed 0
fi

# b. Run S-LDMC on each permuted pair: eQTL sumstats and genotype sample mapping come from the tissue
# (the LD has to match the eQTL sample), while the borzoi effects and borzoi annotations come from the
# randomly drawn partner tissue. Output stems carry both tissues so they cannot collide with the
# matched results in sldmc_results_output_dir.
if false; then
tail -n +2 "$tissue_permuted_pairs_file" | while IFS=$'\t' read -r target_tissue target_sample permuted_tissue permuted_sample; do
	eqtl_sumstats_file=$eqtl_sumstats_dir"eqtl_results_"${target_tissue}"_sumstats.txt.gz"
	genotype_sample_mapping_file=$processed_genotype_data_dir"genotype_sample_mapping_to_"${target_tissue}"_expression_samples.txt"
	borzoi_effect_file=${borzoi_output_dir}${permuted_tissue}"_"${permuted_sample}"_borzoi_effects.txt.gz"
	genotype_stem=$processed_genotype_data_dir"gtex_v9_eqtl_chr"

	annotation_version="default"
	borzoi_annotation_file=${borzoi_output_dir}${permuted_tissue}"_"${permuted_sample}"_annotations_"${annotation_version}".txt.gz"
	sldmc_output_stem=${tissue_permuted_sldmc_results_output_dir}"sldmc_results_eqtl_"${target_tissue}"_borzoi_"${permuted_tissue}"_"${permuted_sample}"_"${annotation_version}
	sbatch run_sldmc.sh $borzoi_effect_file $eqtl_sumstats_file $borzoi_annotation_file $genotype_stem $genotype_sample_mapping_file ${bootstrapped_cross_tissue_gene_sets_dir}"cross_tissue_gene_set_bootstrap_" $sldmc_output_stem $sldmc_code_dir
done
fi






#################
# 7. Visualize results
#################
if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_sldmc_results.R ${sldmc_results_output_dir} $simulation_results_dir $borzoi_gtex_unique_target_names_file $visualize_sldmc_results_dir $annotation_name_file $simulation_oracle_results_dir
fi