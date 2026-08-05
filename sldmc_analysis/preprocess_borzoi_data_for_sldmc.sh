#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 





gtex_v10_pc_genes_gtf="${1}"
borzoi_results_dir="${2}"
borzoi_target_index="${3}"
target_tissue="${4}"
target_sample="${5}"
borzoi_output_dir="${6}"

source ~/.bashrc
conda activate plink_env


echo $target_tissue"_"${target_sample}

python preprocess_borzoi_data_for_sldmc.py $gtex_v10_pc_genes_gtf $borzoi_results_dir $borzoi_target_index $target_tissue $target_sample $borzoi_output_dir
