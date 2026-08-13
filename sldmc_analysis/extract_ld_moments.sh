#!/bin/bash
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=53GB 



borzoi_effect_file="${1}"
eqtl_effects_file="${2}"
borzoi_annotation_file="${3}"
genotype_stem="${4}"
genotype_sample_mapping_file="${5}"
ld_corr_output_stem="${6}"


source ~/.bashrc
conda activate plink_env

python extract_ld_moments.py \
	--est-borzoi-effect-size-file $borzoi_effect_file \
	--est-eqtl-effect-size-file $eqtl_effects_file \
	--variant-gene-annotation-file $borzoi_annotation_file \
	--genotype-plink-filestem $genotype_stem \
	--genotype-sample-mapping-file $genotype_sample_mapping_file \
	--ld-moment-output-stem $ld_corr_output_stem


