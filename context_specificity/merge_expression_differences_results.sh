#!/bin/bash
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



expression_differences_results_dir="${1}"
anno_method="${2}"

source ~/.bashrc
conda activate plink_env


python merge_expression_differences_results.py \
	--expression-differences-results-dir $expression_differences_results_dir \
	--anno-method $anno_method
