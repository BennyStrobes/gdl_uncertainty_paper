#!/bin/bash
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 






source ~/.bashrc
conda activate plink_env


expression_file="${1}"
genotype_sample_mapping_file="${2}"
plink_stem="${3}"
gene_gtf_file="${4}"
new_eqtl_sumstats_output_file="${5}"

echo $expression_file


python simple_eqtl_analysis.py $expression_file $genotype_sample_mapping_file $plink_stem $gene_gtf_file $new_eqtl_sumstats_output_file