#!/bin/bash
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=35GB



genotype_stem="${1}"
sldmc_results_file="${2}"
expr_differences_output_file="${3}"
anno_method="${4}"
borzoi_effect_file1="${5}"
borzoi_annotation_file1="${6}"
genotype_sample_mapping_file1="${7}"
expr_file1="${8}"
borzoi_effect_file2="${9}"
borzoi_annotation_file2="${10}"
genotype_sample_mapping_file2="${11}"
expr_file2="${12}"

source ~/.bashrc
conda activate plink_env


python predict_expression_differences_across_tissue_pairs_in_shared_individuals.py \
	--genotype-stem $genotype_stem \
	--sldmc-results-file $sldmc_results_file \
	--expr-differences-output-file $expr_differences_output_file \
	--anno-method $anno_method \
	--borzoi-effect-file1 $borzoi_effect_file1 \
	--borzoi-annotation-file1 $borzoi_annotation_file1 \
	--genotype-sample-mapping-file1 $genotype_sample_mapping_file1 \
	--expr-file1 $expr_file1 \
	--borzoi-effect-file2 $borzoi_effect_file2 \
	--borzoi-annotation-file2 $borzoi_annotation_file2 \
	--genotype-sample-mapping-file2 $genotype_sample_mapping_file2 \
	--expr-file2 $expr_file2
