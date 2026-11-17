#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 



input_plink_genotype_file_stem="${1}"
genotype_output_dir="${2}"


source ~/.bashrc
conda activate plink_env

mkdir -p ${genotype_output_dir}

# Split into one plink (bed/bim/fam) file per chromosome. No variant filtering.
# Input variant IDs are already GTEx format (chr{chrom}_{pos}_{ref}_{alt}_b38) with A1=ALT / A2=REF, so they are kept as-is.
for chrom_num in {1..22}; do
	echo $chrom_num
	plink2 \
	  --bfile ${input_plink_genotype_file_stem} \
	  --chr ${chrom_num} \
	  --make-bed \
	  --out ${genotype_output_dir%/}/gtex_v8_eqtl_chr${chrom_num}
done


