#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB





gtex_v10_pc_genes_gtf="${1}"
borzoi_results_dir="${2}"
orig_target_index="${3}"
gtex_only_target_index="${4}"
borzoi_gtex_only_target_names_file="${5}"
target_tissue="${6}"
target_sample="${7}"
borzoi_output_dir="${8}"

source ~/.bashrc
conda activate plink_env


echo $target_tissue"_"${target_sample}

python preprocess_rare_variant_borzoi_effects.py $gtex_v10_pc_genes_gtf $borzoi_results_dir $orig_target_index $gtex_only_target_index $borzoi_gtex_only_target_names_file $target_tissue $target_sample $borzoi_output_dir
