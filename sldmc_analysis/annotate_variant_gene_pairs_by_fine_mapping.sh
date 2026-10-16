#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB



borzoi_effect_file="${1}"
fine_mapping_summary_file="${2}"
tissue_name="${3}"
borzoi_annotation_file="${4}"

source ~/.bashrc
conda activate plink_env


echo $borzoi_annotation_file

python annotate_variant_gene_pairs_by_fine_mapping.py \
	--borzoi-effect-file $borzoi_effect_file \
	--fine-mapping-summary-file $fine_mapping_summary_file \
	--tissue-name $tissue_name \
	--borzoi-annotation-file $borzoi_annotation_file \
	--fine-mapping-method "SUSIE" \
	--pip-threshold 0.9
