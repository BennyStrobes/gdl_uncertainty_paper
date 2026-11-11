#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB



sample_size="${1}"
simulation_dir="${2}"


source ~/.bashrc
conda activate plink_env


echo ${sample_size}

python run_no_ld_simulation.py $sample_size $simulation_dir
