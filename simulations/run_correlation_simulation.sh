#!/bin/bash
#SBATCH -t 0-26:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=40GB  


simulation_iter="${1}"
gene_ld_summary_file="${2}"
causal_effect_dir="${3}"
est_eqtl_effect_size_dir="${4}"
est_borzoi_effect_size_dir="${5}"
onek_genomes_plink_filestem="${6}"
inf_output_dir="${7}"
sldmc_code_dir="${8}"

echo "Simulation "${simulation_iter}
source ~/.bashrc
conda activate plink_env


date
####################################################
# Part 1: Simulate causal variant-gene effect sizes
####################################################
echo "PART 1"
causal_variant_gene_effect_size_file=${causal_effect_dir}"sim"${simulation_iter}"_sim_causal_variant_gene_effects.txt.gz"
python simulate_causal_variant_gene_effect_size.py $simulation_iter $gene_ld_summary_file $causal_variant_gene_effect_size_file



####################################################
# Part 2: Simulate estimated (standardized) borzoi sizes
####################################################
echo "PART 2"
n_anno="6"
est_borzoi_standardized_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_est_borzoi_standardized_effects_"${n_anno}"_anno.txt.gz"
sim_variant_gene_annotation_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_annotations.txt.gz"
# SLDMC derives the category file name from the annotation file name, so these two must stay paired
sldmc_variant_gene_annotation_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_annotations_sldmc.txt.gz"
sldmc_annotation_category_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_annotations_sldmc_categories.txt"
python simulate_est_borzoi_effects_for_correlation_experiment.py $causal_variant_gene_effect_size_file $est_borzoi_standardized_effect_size_file ${simulation_iter} $sim_variant_gene_annotation_file $n_anno $sldmc_variant_gene_annotation_file $sldmc_annotation_category_file


####################################################
# Part 2.5: Generate true simulated calibration effect sizes + correlation
####################################################
echo "PART 2.5"
simulation_parameter_summary_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_sim_variant_gene_annotations_"${n_anno}"_true_sim_effect_summary.txt"
python calculate_true_simulated_calibration_effect_sizes_and_correlation.py $est_borzoi_standardized_effect_size_file $sim_variant_gene_annotation_file $causal_variant_gene_effect_size_file $simulation_parameter_summary_file


####################################################
# Part 3: Simulate estimated eqtl effect sizes
####################################################
echo "PART 3"
eqtl_sample_size="489"
est_eqtl_effect_size_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_est_eqtl_effects.txt.gz"
ind_expr_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_individual_expression.txt.gz"
susie_fine_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_susie_fine_mapping.txt.gz"
genotype_sample_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_genotype_sample_mapping.txt"
source ~/.bashrc
conda activate susie
python simulate_eqtl_analysis.py $causal_variant_gene_effect_size_file $est_eqtl_effect_size_file $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $simulation_iter $ind_expr_file $susie_fine_mapping_file $genotype_sample_mapping_file



####################################################
# Part 3.5: Generate true simulated calibration effect sizes + correlation stratified by (annotation, fm_status)
# fm_status: 0 = PIP < thresh (or missing from susie file), 1 = PIP >= thresh
####################################################
echo "PART 3.5"
fm_pip_thresh="0.9"
fm_stratified_simulation_parameter_summary_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_true_fm_stratified_effect_summary.txt"
source ~/.bashrc
conda activate plink_env
python calculate_true_fm_stratified_calibration_effect_sizes_and_correlation.py $est_borzoi_standardized_effect_size_file $sim_variant_gene_annotation_file $causal_variant_gene_effect_size_file $susie_fine_mapping_file $fm_pip_thresh $fm_stratified_simulation_parameter_summary_file


####################################################
# Part 4: Simulate estimated eqtl effect sizes
####################################################
echo "PART 4"
eqtl_sample_size="489"
est_borzoi_effect_size_file=${est_borzoi_effect_size_dir}"sim"${simulation_iter}"_est_borzoi_effects_"${n_anno}"_anno_eqtl_ss_"${eqtl_sample_size}".txt.gz"
source ~/.bashrc
conda activate plink_env
python convert_borzoi_standardized_effects_to_per_allele_effects.py $est_eqtl_effect_size_file $est_borzoi_standardized_effect_size_file $est_borzoi_effect_size_file



####################################################
# Part 5: Run LD corr inference
####################################################
echo "PART 5"
# Updated code
source ~/.bashrc
conda activate sldmc
ld_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_ld_corr_results"
python ${sldmc_code_dir}sldmc.py \
    --est-borzoi-effect-size-file $est_borzoi_effect_size_file \
    --est-eqtl-effect-size-file $est_eqtl_effect_size_file \
    --sim-variant-gene-annotation-file $sldmc_variant_gene_annotation_file \
    --genotype-plink-filestem $onek_genomes_plink_filestem \
    --genotype-sample-mapping-file $genotype_sample_mapping_file \
    --ld-corr-output-stem $ld_corr_output_stem 


####################################################
# Part 5.5: Run LD corr inference with confidently-fine-mapped annotations
# fm_status: 0 = PIP < 0.9 (or missing from susie file), 1 = PIP >= 0.9
# Also includes each simulated annotation crossed with fm_status
####################################################
echo "PART 5.5"
fm_pip_thresh="0.9"
# SLDMC derives the category file name from the annotation file name, so these two must stay paired
fm_status_sldmc_annotation_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_fm_status_annotations_sldmc.txt.gz"
fm_status_sldmc_annotation_category_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_fm_status_annotations_sldmc_categories.txt"
fm_status_ld_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_fm_status_ld_corr_results"
python generate_fine_mapped_sldmc_annotation.py $sldmc_variant_gene_annotation_file $susie_fine_mapping_file $fm_pip_thresh $fm_status_sldmc_annotation_file $fm_status_sldmc_annotation_category_file
python ${sldmc_code_dir}sldmc.py \
    --est-borzoi-effect-size-file $est_borzoi_effect_size_file \
    --est-eqtl-effect-size-file $est_eqtl_effect_size_file \
    --sim-variant-gene-annotation-file $fm_status_sldmc_annotation_file \
    --genotype-plink-filestem $onek_genomes_plink_filestem \
    --genotype-sample-mapping-file $genotype_sample_mapping_file \
    --ld-corr-output-stem $fm_status_ld_corr_output_stem 


####################################################
# Part 6: Run correlations (borzoi vs PMCES) stratified by annotation, in:
# A. Confidently fine-mapped snps (PIP >= 0.9)
# B. All other snps (PIP < 0.9)
# C. All snps regardless of PIP
####################################################
source ~/.bashrc
conda activate plink_env
echo "PART 6"
fm_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_fm_corr_results"
non_fm_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_non_fm_corr_results"
all_snp_pmces_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_all_snp_pmces_corr_results"
python run_fine_map_corr.py $est_borzoi_standardized_effect_size_file $susie_fine_mapping_file $sim_variant_gene_annotation_file $onek_genomes_plink_filestem $fm_corr_output_stem $non_fm_corr_output_stem $all_snp_pmces_corr_output_stem


date

