#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



borzoi_effect_file="${1}"
annotation_name_file="${2}"
borzoi_annotation_filestem="${3}"
eqtl_sumstats_file="${4}"
baselineLD_anno_dir="${5}"
gene_set_anno_file="${6}"

source ~/.bashrc
conda activate plink_env


echo $borzoi_annotation_filestem


##############################
# Version without annotations stratified by magnitude
##############################
magnitude_stratification="False"
borzoi_annotation_file=${borzoi_annotation_filestem}"_default.txt.gz"
python annotate_variant_gene_pairs.py \
	--borzoi-effect-file $borzoi_effect_file \
	--annotation-name-file $annotation_name_file \
	--borzoi-annotation-file $borzoi_annotation_file \
	--eqtl-sumstats-file $eqtl_sumstats_file \
	--baselineLD-anno-dir $baselineLD_anno_dir \
	--gene-set-anno-file $gene_set_anno_file \
	--stratify-by-borzoi-magnitude $magnitude_stratification

##############################
# Version with annotations stratified by magnitude
##############################
magnitude_stratification="True"
borzoi_annotation_file=${borzoi_annotation_filestem}"_magnitude_stratified.txt.gz"
python annotate_variant_gene_pairs.py \
	--borzoi-effect-file $borzoi_effect_file \
	--annotation-name-file $annotation_name_file \
	--borzoi-annotation-file $borzoi_annotation_file \
	--eqtl-sumstats-file $eqtl_sumstats_file \
	--baselineLD-anno-dir $baselineLD_anno_dir \
	--gene-set-anno-file $gene_set_anno_file \
	--stratify-by-borzoi-magnitude $magnitude_stratification
