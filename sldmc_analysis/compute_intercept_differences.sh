#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



ld_corr_output_file_list="${1}"
meta_analyzed_output_stem="${2}"
annotation_version="${3}"


source ~/.bashrc
conda activate plink_env

date

python compute_intercept_differences.py \
	--ld-corr-output-file-list $ld_corr_output_file_list \
	--meta-analyzed-output-stem $meta_analyzed_output_stem \
	--annotation-version $annotation_version

date
