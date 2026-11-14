#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



sldmc_summary_file="${1}"
borzoi_results_file="${2}"
afr_expression_file="${3}"
plink_genotype_stem="${4}"
afr_genotype_sample_mapping_file="${5}"
eur_genotype_sample_mapping_file="${6}"
expression_correlation_output_file="${7}"


source ~/.bashrc
conda activate plink_env


echo ${expression_correlation_output_file}

python african_personalized_expression_correlations_per_tissue.py \
	--sldmc-summary-file $sldmc_summary_file \
	--borzoi-results-file $borzoi_results_file \
	--expression-file $afr_expression_file \
	--plink-genotype-stem $plink_genotype_stem \
	--genotype-sample-mapping-file $afr_genotype_sample_mapping_file \
	--eur-genotype-sample-mapping-file $eur_genotype_sample_mapping_file \
	--expression-correlation-output-file $expression_correlation_output_file
