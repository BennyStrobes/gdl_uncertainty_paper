#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



sldmc_summary_file="${1}"
rare_variant_borzoi_results_file="${2}"
common_variant_borzoi_results_file="${3}"
expression_file="${4}"
plink_genotype_stem="${5}"
genotype_sample_mapping_file="${6}"
rare_variant_carrier_stem="${7}"
eur_vcf_sample_list="${8}"
maf_threshold="${9}"
expression_correlation_output_file="${10}"


source ~/.bashrc
conda activate plink_env


echo ${expression_correlation_output_file}

python rare_variant_personalized_expression_correlations_per_tissue.py \
	--sldmc-summary-file $sldmc_summary_file \
	--rare-variant-borzoi-results-file $rare_variant_borzoi_results_file \
	--common-variant-borzoi-results-file $common_variant_borzoi_results_file \
	--expression-file $expression_file \
	--plink-genotype-stem $plink_genotype_stem \
	--genotype-sample-mapping-file $genotype_sample_mapping_file \
	--rare-variant-carrier-stem $rare_variant_carrier_stem \
	--eur-vcf-sample-list $eur_vcf_sample_list \
	--maf-threshold $maf_threshold \
	--expression-correlation-output-file $expression_correlation_output_file
