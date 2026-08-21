#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 





residual_expression_file="${1}"
plink_fam_file="${2}"
genotype_mapping_file="${3}"



source ~/.bashrc
conda activate plink_env


python create_genotype_sample_mapping.py $residual_expression_file $plink_fam_file $genotype_mapping_file
