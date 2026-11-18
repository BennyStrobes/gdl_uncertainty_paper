#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB



input_expression_file="${1}"
input_covariate_file="${2}"
afr_residual_expression_file="${3}"
eur_residual_expression_file="${4}"
genes_gtf="${5}"
gtex_donor_info_file="${6}"
genotype_fam_file="${7}"
borzoi_results_dir="${8}"


source ~/.bashrc
conda activate plink_env



python residualize_expression_all_sample_fit.py $input_expression_file $input_covariate_file $afr_residual_expression_file $eur_residual_expression_file $genes_gtf $gtex_donor_info_file $genotype_fam_file $borzoi_results_dir
