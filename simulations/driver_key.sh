#!/bin/bash
#SBATCH -t 0-3:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                          # Partition to run in
#SBATCH --mem=10GB  

######################
# Input data
######################
# File stem to 1 thousand genomes plink files
onek_genomes_plink_filestem="/lab-share/CHIP-Strober-e2/Public/1000G_Phase3/hg38/1000G.EUR.hg38."

# Gene annotation file
gene_annotation_file="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"

# Code directory for SLDMC
sldmc_code_dir="/lab-share/CHIP-Strober-e2/Public/ben/SLDMC/"

######################
# Output data
######################
# Output root
output_root="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/simulations/"

# Gene info directory
gene_info_dir=${output_root}"gene_info/"

# simulated correlation causal eqtl effects directory
correlation_causal_effect_dir=${output_root}"corr_sim_causal_eqtl_effects/"

correlation_borzoi_est_effect_dir=${output_root}"corr_sim_borzoi_eqtl_effects/"

correlation_est_eqtl_effects_dir=${output_root}"corr_sim_est_eqtl_effect_size_dir/"

correlation_inference_results_dir=${output_root}"corr_sim_inference_results/"

correlation_visualization_dir=${output_root}"corr_visualization/"




######################
# Run analysis
######################

# First generate gene info including list of genes and variant-to-gene links
gene_summary_file=$gene_info_dir"gene_summary.txt"
if false; then
source ~/.bashrc
conda activate plink_env
python generate_simulation_gene_info.py $gene_annotation_file $onek_genomes_plink_filestem $gene_info_dir
fi


######################
# Run simulations based to get correlation of predicted and observed causal effects
######################
if false; then
for simulation_iter in {1..50}
do
	sbatch run_correlation_simulation.sh $simulation_iter $gene_summary_file $correlation_causal_effect_dir $correlation_est_eqtl_effects_dir $correlation_borzoi_est_effect_dir $onek_genomes_plink_filestem $correlation_inference_results_dir $sldmc_code_dir
done
fi


######################
# Iterate simulations across multiple eQTL sample sizes
######################
if false; then
for simulation_iter in {1..50}
do
	sbatch run_multi_sample_size_correlation_simulation.sh $simulation_iter $gene_summary_file $correlation_causal_effect_dir $correlation_est_eqtl_effects_dir $correlation_borzoi_est_effect_dir $onek_genomes_plink_filestem $correlation_inference_results_dir $sldmc_code_dir
done
fi




if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_corr_simulation_results.R $correlation_inference_results_dir $correlation_borzoi_est_effect_dir $correlation_visualization_dir
fi