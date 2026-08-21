#!/bin/bash
#SBATCH -t 0-0:12                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 



input_expression_file="${1}"
input_covariate_file="${2}"
residual_expression_file="${3}"
gtex_v10_pc_genes_gtf="${4}"
gtex_donor_info_file="${5}"
genotype_fam_file="${6}"


source ~/.bashrc
conda activate plink_env



python residualize_expression.py $input_expression_file $input_covariate_file $residual_expression_file $gtex_v10_pc_genes_gtf $gtex_donor_info_file $genotype_fam_file