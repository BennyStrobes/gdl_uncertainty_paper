#!/bin/bash
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=53GB 



borzoi_effect_file="${1}"
eqtl_effects_file="${2}"
borzoi_annotation_file="${3}"
genotype_stem="${4}"
genotype_sample_mapping_file="${5}"
bootstrapped_cross_tissue_gene_sets_file_stem="${6}"
ld_corr_output_stem="${7}"
sldmc_code_dir="${8}"


source ~/.bashrc
conda activate plink_env

date
echo $ld_corr_output_stem

python ${sldmc_code_dir}sldmc.py \
	--est-borzoi-effect-size-file $borzoi_effect_file \
	--est-eqtl-effect-size-file $eqtl_effects_file \
	--sim-variant-gene-annotation-file $borzoi_annotation_file \
	--genotype-plink-filestem $genotype_stem \
	--genotype-sample-mapping-file $genotype_sample_mapping_file \
	--ld-corr-output-stem $ld_corr_output_stem \
	--bootstrapped-gene-set-filestem ${bootstrapped_cross_tissue_gene_sets_file_stem} \
    --weighted "True"


date
