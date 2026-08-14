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
# Assumes Part 1 and 2 have already been run and files exist
####################################################



####################################################
# Part 3: Simulate estimated eqtl effect sizes
####################################################
echo "PART 3"
eqtl_sample_size="450"
est_eqtl_effect_size_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_est_eqtl_effects.txt.gz"
ind_expr_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_individual_expression.txt.gz"
susie_fine_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_susie_fine_mapping.txt.gz"
genotype_sample_mapping_file=${est_eqtl_effect_size_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_genotype_sample_mapping.txt"
source ~/.bashrc
conda activate susie
python simulate_eqtl_analysis.py $causal_variant_gene_effect_size_file $est_eqtl_effect_size_file $gene_ld_summary_file $onek_genomes_plink_filestem $eqtl_sample_size $simulation_iter $ind_expr_file $susie_fine_mapping_file $genotype_sample_mapping_file


####################################################
# Part 4: Simulate estimated eqtl effect sizes
####################################################
echo "PART 4"
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
# Part 6: Run correlations based on only fine-mapped snps
####################################################
source ~/.bashrc
conda activate plink_env
echo "PART 6"
fm_corr_output_stem=${inf_output_dir}"sim"${simulation_iter}"_sim_eqtl_ss_"${eqtl_sample_size}"_"${n_anno}"_anno_fm_corr_results"
python run_fine_map_corr.py $est_borzoi_standardized_effect_size_file $susie_fine_mapping_file $sim_variant_gene_annotation_file $onek_genomes_plink_filestem $fm_corr_output_stem


date

