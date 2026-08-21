#################
# Input data
#################

# Directory containing results of borzoi runs
borzoi_results_dir="/lab-share/CHIP-Strober-e2/Public/ben/borzoi_genome_wide_run/genome_wide/borzoi_predictions/"

# Directory containing borzoi gtex target indices and names
borzoi_gtex_target_names_file=${borzoi_results_dir}"targets_gtex_only_ordered.txt"
borzoi_gtex_independent_target_names_file=${borzoi_results_dir}"targets_gtex_only_independent_ordered.txt"
borzoi_gtex_independent_target_names_file=${borzoi_results_dir}"targets_gtex_only_unique_ordered.txt"



# Directory containing genotype data
plink2_genotype_data_dir="/lab-share/CHIP-Strober-e2/Public/ben/process_gtex_genotype_data/processed_genotype/"

# Gtex donor attributes (has ancestry)
gtex_donor_info_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_899Indiv_QC_TABLE.tsv"

# GTex v8 per tissue expression
gtex_v8_expression_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_v8_eQTL_expression_matrices/"

# GTEx v8 covariate dir
gtex_v8_covariate_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_v8_eQTL_covariates/"



# Directory containing eQTL summary statistics
eqtl_sumstats_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/eqtl_sumstats/"

# Gtex v10 protein coding genes
gtex_v10_pc_genes_gtf="/lab-share/CHIP-Strober-e2/Public/gene_annotation_files/gencode.v39.gtex.protein_coding.genes.gtf"


#################
# Output data
#################
output_root_dir="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/gtex_eqtl_expression_processing/"

# Directory containing residualized expression
residualized_expression_dir=${output_root_dir}"residualized_expression/"

# Directory containing processed genotype
plink_genotype_data_dir=${output_root_dir}"plink_processed_genotype/"

# Directory continaing processed eQTL results
eQTL_results_dir=${output_root_dir}"eqtl_results/"


#################
# Run code
#################



############################
# 0. Convert plink2 to plink genotype
############################
if false; then
sbatch convert_plink2_to_plink_genotype_files.sh $plink2_genotype_data_dir $plink_genotype_data_dir
fi


############################
# 1. Residualize expression (ie. remove covariates)
############################
genotype_fam_file=$plink_genotype_data_dir"gtex_v9_eqtl_chr1.fam"
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do

	echo $gtex_tissue
	input_expression_file=${gtex_v8_expression_dir}${gtex_tissue}".v8.normalized_expression.bed.gz"
	input_covariate_file=${gtex_v8_covariate_dir}${gtex_tissue}".v8.covariates.txt"

	residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8.residualized_expression_renormalized.bed"
	sbatch residualize_expression.sh $input_expression_file $input_covariate_file $residual_expression_file $gtex_v10_pc_genes_gtf $gtex_donor_info_file $genotype_fam_file
done
fi


############################
# 2. Create mapping of ordered genotype indices for each tissue
############################
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	echo $gtex_tissue
	residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8.residualized_expression_renormalized.bed"
	genotype_mapping_file=$plink_genotype_data_dir"genotype_sample_mapping_to_"${gtex_tissue}"_expression_samples.txt"
	sh create_genotype_sample_mapping.sh $residual_expression_file $plink_genotype_data_dir"gtex_v9_eqtl_chr1.fam" $genotype_mapping_file
done
fi






############################
# 3. Run Simple eQTL analsysis
############################
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8.residualized_expression_renormalized.bed"
	genotype_mapping_file=$plink_genotype_data_dir"genotype_sample_mapping_to_"${gtex_tissue}"_expression_samples.txt"

	eqtl_results_output_file=$eQTL_results_dir"eqtl_results_"${gtex_tissue}"_sumstats.txt.gz"
	sbatch simple_eqtl_analysis.sh $residual_expression_file $genotype_mapping_file $plink_genotype_data_dir"gtex_v9_eqtl_chr" $eqtl_sumstats_dir${gtex_tissue}".v10.allpairs.chr" $eqtl_results_output_file $gtex_v10_pc_genes_gtf
done
fi






