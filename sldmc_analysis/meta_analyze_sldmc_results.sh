#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



ld_corr_output_file_list="${1}"
meta_analyzed_ld_corr_output_stem="${2}"


source ~/.bashrc
conda activate plink_env

date

python meta_analyze_sldmc_results.py \
	--ld-corr-output-file-list $ld_corr_output_file_list \
	--meta-analyzed-ld-corr-output-stem $meta_analyzed_ld_corr_output_stem \
	--correlation-meta-method "recompute_from_components"


if false; then
# Alt approach that averages correlations across tissues instead of recomputing from components. This is not the approach we used in the paper, but it is included here for completeness.
python meta_analyze_sldmc_results.py \
	--ld-corr-output-file-list $ld_corr_output_file_list \
	--meta-analyzed-ld-corr-output-stem $meta_analyzed_ld_corr_output_stem"_avg_correlation" \
	--correlation-meta-method "average_correlation"
fi

date
