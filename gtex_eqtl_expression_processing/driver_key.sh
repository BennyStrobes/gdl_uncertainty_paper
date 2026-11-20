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
# Suffix is .bim, .bed, .fam
input_plink_genotype_file_stem="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/josh_popp_processed_eqtl_expression_data/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01"

# Gtex donor attributes (has ancestry)
gtex_donor_info_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_899Indiv_QC_TABLE.tsv"

# GTex v8 per tissue expression
gtex_v8_expression_dir="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/josh_popp_processed_eqtl_expression_data/expression/"

# GTEx v8 covariate dir
gtex_v8_covariate_dir="/lab-share/CHIP-Strober-e2/Public/ben/gdl_uncertainty_paper/josh_popp_processed_eqtl_expression_data/covariates/"

# GTF file used by Borzoi
gtf_file="/home/ch271704/tools/borzoi/examples/hg38/genes/gencode41/gencode41_basic_nort.gtf"


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
# 0. Convert plink genome wide genotype data into one plink file per chromosome
############################
if false; then
sh split_genotype_data_into_chromosomes.sh $input_plink_genotype_file_stem $plink_genotype_data_dir
fi



###################
# Part A (main analysis): Generate data restricted to European ancestry individuals
###################

############################
# 1. Residualize expression (ie. remove covariates)
############################
genotype_fam_file=$plink_genotype_data_dir"gtex_v8_eqtl_chr1.fam"
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do

	echo $gtex_tissue
	input_expression_file=${gtex_v8_expression_dir}${gtex_tissue}".v8.normalized_expression.bed.gz"
	input_covariate_file=${gtex_v8_covariate_dir}${gtex_tissue}".v8.covariates.txt"

	residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8.residualized_expression_renormalized.bed"
	sbatch residualize_expression.sh $input_expression_file $input_covariate_file $residual_expression_file $gtf_file $gtex_donor_info_file $genotype_fam_file "EUR" $borzoi_results_dir
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
	sh create_genotype_sample_mapping.sh $residual_expression_file $plink_genotype_data_dir"gtex_v8_eqtl_chr1.fam" $genotype_mapping_file
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
	sbatch simple_eqtl_analysis.sh $residual_expression_file $genotype_mapping_file $plink_genotype_data_dir"gtex_v8_eqtl_chr" $gtf_file $eqtl_results_output_file
done
fi





###################
# Part B (Secondary Analysis): Generate data restricted to African ancestry individuals
###################

############################
# 1. Residualize expression (ie. remove covariates)
# Covariate effects (full GTEx covariate set) are fit on ALL samples of the tissue, then residuals are
# subset to AFR (and, separately, EUR) individuals and standardized within ancestry.
# The EUR output (.v8_eur_allfit.) is a separate file; Part A's EUR files are not touched.
############################
genotype_fam_file=$plink_genotype_data_dir"gtex_v8_eqtl_chr1.fam"
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do

	echo $gtex_tissue
	input_expression_file=${gtex_v8_expression_dir}${gtex_tissue}".v8.normalized_expression.bed.gz"
	input_covariate_file=${gtex_v8_covariate_dir}${gtex_tissue}".v8.covariates.txt"

	afr_residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8_afr.residualized_expression_renormalized.bed"
	eur_allfit_residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8_eur_allfit.residualized_expression_renormalized.bed"
	sbatch residualize_expression_all_sample_fit.sh $input_expression_file $input_covariate_file $afr_residual_expression_file $eur_allfit_residual_expression_file $gtf_file $gtex_donor_info_file $genotype_fam_file $borzoi_results_dir
done
fi



############################
# 2. Create mapping of ordered genotype indices for each tissue
############################
if false; then
tail -n +2 "$borzoi_gtex_independent_target_names_file" | while IFS=$'\t' read -r orig_target_index borzoi_target_index target_identifier target_description gtex_tissue; do
	echo $gtex_tissue
	residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8_afr.residualized_expression_renormalized.bed"
	genotype_mapping_file=$plink_genotype_data_dir"genotype_sample_mapping_to_"${gtex_tissue}"_afr_expression_samples.txt"
	sh create_genotype_sample_mapping.sh $residual_expression_file $plink_genotype_data_dir"gtex_v8_eqtl_chr1.fam" $genotype_mapping_file

	# EUR (all-sample covariate fit): should have identical samples/order to Part A EUR, so Part A's EUR mapping can be reused. Verify.
	eur_allfit_residual_expression_file=${residualized_expression_dir}${gtex_tissue}".v8_eur_allfit.residualized_expression_renormalized.bed"
	eur_allfit_genotype_mapping_file=$plink_genotype_data_dir"genotype_sample_mapping_to_"${gtex_tissue}"_eur_allfit_expression_samples.txt"
	sh create_genotype_sample_mapping.sh $eur_allfit_residual_expression_file $plink_genotype_data_dir"gtex_v8_eqtl_chr1.fam" $eur_allfit_genotype_mapping_file
	eur_genotype_mapping_file=$plink_genotype_data_dir"genotype_sample_mapping_to_"${gtex_tissue}"_expression_samples.txt"
	if cmp -s $eur_allfit_genotype_mapping_file $eur_genotype_mapping_file; then
		echo "EUR allfit sample mapping matches Part A"
	else
		echo "WARNING: EUR allfit sample mapping differs from Part A for "$gtex_tissue
	fi
done
fi


