###################
# Input data
###################






#####################
# Output data
#####################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/cross_context_heritability/"

simulation_dir="${output_root}no_ld_simulation/"

simulation_visualization_dir="${output_root}simulation_visualization/"


#######################
# Scripts
########################


#####################
# No LD simulation
sample_size="600"
if false; then
sh run_no_ld_simulation.sh $sample_size $simulation_dir
fi

#####################
# Visualize No LD simulation
if false; then
source ~/.bashrc
conda activate plink_env
Rscript visualize_no_ld_simulation_results.R $simulation_dir $sample_size $simulation_visualization_dir
fi
