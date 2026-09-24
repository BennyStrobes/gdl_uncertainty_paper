#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



sldmc_summary_file="${1}"
borzoi_results_file="${2}"
expression_file="${3}"
plink_genotype_stem="${4}"
genotype_sample_mapping_file="${5}"
expression_correlation_output_file="${6}"


source ~/.bashrc
conda activate plink_env


echo ${expression_correlation_output_file}

python personalized_expression_correlations_per_tissue.py \
	--sldmc-summary-file $sldmc_summary_file \
	--borzoi-results-file $borzoi_results_file \
	--expression-file $expression_file \
	--plink-genotype-stem $plink_genotype_stem \
	--genotype-sample-mapping-file $genotype_sample_mapping_file \
	--expression-correlation-output-file $expression_correlation_output_file
