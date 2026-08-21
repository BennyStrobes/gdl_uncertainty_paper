#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=10GB 



plink2_genotype_data_dir="${1}"
genotype_output_dir="${2}"


source ~/.bashrc
conda activate plink_env



for chrom_num in {1..22}; do
	echo $chrom_num
	plink2 \
	  --pfile ${plink2_genotype_data_dir%/}/gtex_v9_eqtl_chr${chrom_num} \
	  --set-all-var-ids chr@_#_\$r_\$a_b38 \
	  --make-bed \
	  --out ${genotype_output_dir%/}/gtex_v9_eqtl_chr${chrom_num}
done