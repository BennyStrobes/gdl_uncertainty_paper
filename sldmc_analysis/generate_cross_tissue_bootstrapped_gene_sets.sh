#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB



eqtl_sumstats_dir="${1}"
borzoi_output_dir="${2}"
borzoi_gtex_independent_target_names_file="${3}"
bootstrapped_cross_tissue_gene_sets_dir="${4}"


source ~/.bashrc
conda activate plink_env


python generate_cross_tissue_bootstrapped_gene_sets.py $eqtl_sumstats_dir $borzoi_output_dir $borzoi_gtex_independent_target_names_file $bootstrapped_cross_tissue_gene_sets_dir
