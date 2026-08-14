import numpy as np
import sys
import pdb
import gzip
import argparse
from pandas_plink import read_plink


def create_mapping_from_gene_id_to_per_allele_borzoi_effects(borzoi_effect_file):
    """
    Load per-allele Borzoi effect sizes for each variant-gene pair.
    Returns mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, per_allele_effect)
    """
    f = gzip.open(borzoi_effect_file, 'rt')
    mapping = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        gene_id = data[0]
        var_id = data[1]
        chrom_num = data[2]
        snp_pos = data[3]
        a0 = data[4]
        a1 = data[5]
        if a0 == a1:
            print('assumption eroroor')
            pdb.set_trace()
        per_allele_effect = float(data[6])

        if gene_id not in mapping:
            mapping[gene_id] = {}
        if var_id in mapping[gene_id]:
            print('variatn repeat assumption erororo')
            pdb.set_trace()

        mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, per_allele_effect)
    f.close()
    return mapping


def create_mapping_from_variant_id_to_genotype_index(ordered_snps):
    mapping = {}

    n_snps = len(ordered_snps)
    for snp_iter in range(n_snps):
        snp_name = ordered_snps[snp_iter]

        if snp_name in mapping:
            print('asssumption erororo')
            pdb.set_trace()
        mapping[snp_name] = snp_iter

    return mapping


def create_mapping_from_variant_id_to_snp_info(snp_array, a0_arr, a1_arr, chrom_arr, pos_arr):
    if len(snp_array) != len(a0_arr):
        print('assumption eorroro')
        pdb.set_trace()
    if len(snp_array) != len(a1_arr):
        print('assumption eorroro')
        pdb.set_trace()

    dicti = {}

    for ii, snp_id in enumerate(snp_array):
        if snp_id in dicti:
            print('assumpationoenroer')
            pdb.set_trace()
        dicti[snp_id] = (a0_arr[ii], a1_arr[ii], chrom_arr[ii], pos_arr[ii])
    return dicti


def extract_ordered_cis_variants_for_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_borzoi_effects):
    # Ordered list of the gene's cis variants that are present in the genotype data with consistent alleles
    unique_vars = np.unique([*var_to_borzoi_effects])
    final_vars = []
    for var in unique_vars:
        if var not in rsid_to_genotype_index:
            continue
        geno_alleles = (rsid_to_snp_info[var][0], rsid_to_snp_info[var][1])
        borzoi_alleles = var_to_borzoi_effects[var][4:6]
        if set(geno_alleles) != set(borzoi_alleles):
            continue
        final_vars.append(var)
    return np.asarray(final_vars)


def create_mapping_from_gene_id_to_variant_gene_annotations(borzoi_annotation_file, anno_method):
    # Columns 6+ are one integer column per annotation, holding the category index
    # the variant-gene pair falls in (-1 if the pair falls in no category of that annotation).
    # Only the anno_method column is kept in memory.
    f = gzip.open(borzoi_annotation_file, 'rt')
    mapping = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            anno_names = np.asarray(data[6:])
            # Column of the annotation file corresponding to the target annotation
            if np.sum(anno_names == anno_method) != 1:
                print('assumption eroror: ' + anno_method + ' not found in annotation file')
                pdb.set_trace()
            target_anno_column = 6 + np.where(anno_names == anno_method)[0][0]
            continue
        gene_id = data[0]
        var_id = data[1]
        chrom_num = data[2]
        snp_pos = data[3]
        a0 = data[4]
        a1 = data[5]
        if a0 == a1:
            print('assumption eroroor')
            pdb.set_trace()
        category_index = int(data[target_anno_column])
        if gene_id not in mapping:
            mapping[gene_id] = {}
        if var_id in mapping[gene_id]:
            print('variatn repeat assumption erororo')
            pdb.set_trace()
        mapping[gene_id][var_id] = (gene_id, var_id, chrom_num, snp_pos, a0, a1, category_index)
    f.close()
    return mapping, anno_names


def extract_annotation_categories(annotation_category_file, anno_names):
    # Companion file to the annotation file. Columns (tab-separated, with header):
    # anno_name  source  category_index  category_name
    # Returns, per annotation, the ordered list of its category names.
    anno_name_to_category_names = {}
    anno_name_to_source = {}
    f = open(annotation_category_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        anno_name = data[0]
        source = data[1]
        category_index = int(data[2])
        category_name = data[3]
        if anno_name not in anno_name_to_category_names:
            anno_name_to_category_names[anno_name] = []
            anno_name_to_source[anno_name] = source
        # Categories are written in index order, so appending keeps them aligned to the index
        if category_index != len(anno_name_to_category_names[anno_name]):
            print('assumption eroror: categories not in index order for ' + anno_name)
            pdb.set_trace()
        anno_name_to_category_names[anno_name].append(category_name)
    f.close()

    # Every annotation column in the annotation file needs an entry in the category file
    for anno_name in anno_names:
        if anno_name not in anno_name_to_category_names:
            print('assumption eroror: ' + anno_name + ' missing from category file')
            pdb.set_trace()

    return anno_name_to_category_names, anno_name_to_source


def extract_sldmc_calibration_and_uncertainty(sldmc_results_file, anno_method):
    """
    Extract per-category calibration slopes and effect-prediction uncertainty from the SLDMC results file.
    Returns mapping[category_name] = {output_name: point_estimate} for the requested annotation, with
    output_name in (calibration_slope, per_snp_eqtl_h2, borzoi_variance, unstandardized_borzoi_variance, correlation),
    plus residual_variance = per_snp_eqtl_h2 - calibration_slope^2 * borzoi_variance (the variance of
    true standardized eQTL effects around the calibrated Borzoi prediction).
    """
    f = open(sldmc_results_file)
    mapping = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        annotation_name = data[0]
        if annotation_name != anno_method:
            continue
        category_name = data[1]
        output_name = data[2]
        point_estimate = float(data[3])

        if category_name not in mapping:
            mapping[category_name] = {}
        if output_name in mapping[category_name]:
            print('repeated output name assumption eroror')
            pdb.set_trace()
        mapping[category_name][output_name] = point_estimate
    f.close()

    if len(mapping) == 0:
        print('no sldmc results found for annotation ' + anno_method)
        pdb.set_trace()

    for category_name in mapping:
        stats = mapping[category_name]
        stats['residual_variance'] = stats['per_snp_eqtl_h2'] - np.square(stats['calibration_slope'])*stats['borzoi_variance']

    return mapping


def get_observed_and_predicted_expression_for_all_gene_individual_pairs_in_a_tissue(genotype_stem, sldmc_results_file, anno_method, borzoi_effect_file, borzoi_annotation_file, genotype_sample_mapping_file, expr_file):
    """
    For a given tissue, get the observed and predicted expression for all gene-individual pairs.
    """

    # Extract per-category calibration slopes and effect-prediction uncertainty (standardized-effect units)
    category_name_to_sldmc_stats = extract_sldmc_calibration_and_uncertainty(sldmc_results_file, anno_method)

    # Load per variant-gene pair category indices of the target annotation
    gene_id_to_variant_gene_anno, anno_names = create_mapping_from_gene_id_to_variant_gene_annotations(borzoi_annotation_file, anno_method)

    # Load per-allele Borzoi effect sizes for each variant-gene pair
    gene_id_to_per_allele_borzoi_effects = create_mapping_from_gene_id_to_per_allele_borzoi_effects(borzoi_effect_file)

    # Ordered category names of the target annotation (maps category index -> category name)
    annotation_category_file = borzoi_annotation_file.split('.txt.gz')[0] + '_categories.txt'
    anno_name_to_category_names, anno_name_to_source = extract_annotation_categories(annotation_category_file, anno_names)
    ordered_category_names = anno_name_to_category_names[anno_method]

    # Every category of the target annotation needs sldmc calibration/uncertainty stats
    for category_name in ordered_category_names:
        if category_name not in category_name_to_sldmc_stats:
            print('assumption eroror: no sldmc stats for category ' + category_name)
            pdb.set_trace()

    # Genotype-column indices aligning genotype samples to the expression sample ordering
    genotype_sample_indices = (np.loadtxt(genotype_sample_mapping_file)).astype(int)

    # Initialize output objects
    gene_individual_to_results = {}
    ordered_sample_names = None

    # Loop through chromsomes
    for chrom_num in range(1, 23):
        print(chrom_num, flush=True)

        ##################################
        # Load in per-chrom-genotype data
        ##################################
        # string of chromosome name
        chrom_string = 'chr' + str(chrom_num)
        # Load in chromosome plink data
        (bim, fam, G) = read_plink(genotype_stem + str(chrom_num))
        # Create mapping from variant id to index
        rsid_to_genotype_index = create_mapping_from_variant_id_to_genotype_index(np.asarray(bim['snp']))
        # Create mapping from rsid to a0, a1
        rsid_to_snp_info = create_mapping_from_variant_id_to_snp_info(np.asarray(bim['snp']), np.asarray(bim['a0']), np.asarray(bim['a1']), np.asarray(bim['chrom']), np.asarray(bim['pos']))

        ##################################
        # Loop through genes on this chromosome
        # (rows of the expression file)
        ##################################
        f = open(expr_file)
        head_count = 0
        for line in f:
            line = line.rstrip()
            data = line.split('\t')
            if head_count == 0:
                head_count = head_count + 1
                # Sample ids of the expression file (columns 4+, aligned with all per-individual vectors)
                ordered_sample_names = np.asarray(data[4:])
                continue

            # Filter out genes not on this chromosome
            line_chrom = data[0]
            if line_chrom != chrom_string:
                continue

            # Expression file gene ids are versioned; borzoi effect file uses short gene ids
            full_gene_id = data[3]
            short_gene_id = full_gene_id.split('.')[0]
            if short_gene_id not in gene_id_to_per_allele_borzoi_effects:
                continue
            var_to_borzoi_effects = gene_id_to_per_allele_borzoi_effects[short_gene_id]

            # Gene also needs variant-gene annotations
            if short_gene_id not in gene_id_to_variant_gene_anno:
                continue
            var_to_anno = gene_id_to_variant_gene_anno[short_gene_id]

            # Observed expression for the gene (one value per individual)
            expr_vec = np.asarray(data[4:]).astype(float)

            # Extract ordered list of the gene's cis variants
            ordered_cis_variants = extract_ordered_cis_variants_for_gene(rsid_to_genotype_index, rsid_to_snp_info, var_to_borzoi_effects)
            if len(ordered_cis_variants) == 0:
                continue

            # Genotype rows for the cis variants, checking allele orientation against borzoi alleles
            cis_genotype_indices = []
            for cis_variant in ordered_cis_variants:
                cis_genotype_indices.append(rsid_to_genotype_index[cis_variant])
                snp_info = rsid_to_snp_info[cis_variant]
                # After the 2.0 - dosage flip below, dosages count the plink a0 allele, so the
                # borzoi effect allele (a1) must be the plink a0 allele
                geno_allele_alt_string = snp_info[1] + '_' + snp_info[0]
                borzoi_allele_string = var_to_borzoi_effects[cis_variant][4] + '_' + var_to_borzoi_effects[cis_variant][5]
                if geno_allele_alt_string != borzoi_allele_string:
                    print('allele orientation assumption eroror')
                    pdb.set_trace()

            # Extract genotype matrix (variants X individuals), columns aligned with expr_vec
            cis_genotype_indices = np.asarray(cis_genotype_indices)
            geno_mat = (G[cis_genotype_indices, :].compute())[:, genotype_sample_indices]

            # Mean impute missing genotypes
            row_means = np.nanmean(geno_mat, axis=1)
            nan_rows, nan_cols = np.where(np.isnan(geno_mat))
            geno_mat[nan_rows, nan_cols] = row_means[nan_rows]

            # Convert to dosage of the borzoi effect allele
            geno_mat = 2.0 - geno_mat

            if geno_mat.shape[1] != len(expr_vec):
                print('sample alignment assumption eroror')
                pdb.set_trace()

            ##################################
            # Standardized genotype and standardized borzoi effects
            ##################################
            # Genotype standard deviation of each variant
            geno_sdevs = np.std(geno_mat, axis=1)
            # Variants with no genotype variance in this sample cannot be standardized
            degenerate_indices = (geno_sdevs > 0.0) == False

            # Standardized genotype (mean zero, unit variance per variant)
            std_geno_mat = geno_mat - np.mean(geno_mat, axis=1)[:, None]
            std_geno_mat[degenerate_indices == False, :] = std_geno_mat[degenerate_indices == False, :]/(geno_sdevs[degenerate_indices == False][:, None])

            # Scale per-allele borzoi effects into standardized-genotype units (the units the SLDMC moments were computed in)
            per_allele_borzoi_effect_vec = []
            for cis_variant in ordered_cis_variants:
                per_allele_borzoi_effect_vec.append(var_to_borzoi_effects[cis_variant][6])
            per_allele_borzoi_effect_vec = np.asarray(per_allele_borzoi_effect_vec)
            std_borzoi_effect_vec = per_allele_borzoi_effect_vec*geno_sdevs
            std_borzoi_effect_vec[degenerate_indices] = np.nan

            ##################################
            # Per-variant calibration factor and prediction uncertainty
            # (based on the variant-gene pair's annotation category)
            ##################################
            calibration_slope_vec = []
            residual_variance_vec = []
            for cis_variant in ordered_cis_variants:
                if cis_variant not in var_to_anno:
                    calibration_slope_vec.append(np.nan)
                    residual_variance_vec.append(np.nan)
                    continue
                if var_to_anno[cis_variant][4] != var_to_borzoi_effects[cis_variant][4]:
                    print('annotation alllele assumption erororo')
                    pdb.set_trace()
                category_index = var_to_anno[cis_variant][6]
                if category_index == -1:
                    # Variant-gene pair falls in no category of the annotation -> no calibration info
                    calibration_slope_vec.append(np.nan)
                    residual_variance_vec.append(np.nan)
                    continue
                category_name = ordered_category_names[category_index]
                sldmc_stats = category_name_to_sldmc_stats[category_name]
                calibration_slope_vec.append(sldmc_stats['calibration_slope'])
                residual_variance_vec.append(sldmc_stats['residual_variance'])
            calibration_slope_vec = np.asarray(calibration_slope_vec)
            residual_variance_vec = np.asarray(residual_variance_vec)

            ##################################
            # Per-individual observed expression, predicted expression, and prediction uncertainty
            ##################################
            # Variants with no calibration info or no standardized effect contribute nothing to the prediction
            valid_variant_indices = (np.isnan(std_borzoi_effect_vec) | np.isnan(calibration_slope_vec)) == False

            # Calibrated standardized effect of each valid variant
            calibrated_std_borzoi_effect_vec = calibration_slope_vec[valid_variant_indices]*std_borzoi_effect_vec[valid_variant_indices]

            # Predicted expression of each individual
            predicted_expr_vec = np.dot(np.transpose(std_geno_mat[valid_variant_indices, :]), calibrated_std_borzoi_effect_vec)

            # Uncertainty variance in predicted expression of each individual
            # (assumes independent effect residuals across variants: sum_v std_geno_iv^2 * residual_variance_v)
            predicted_expr_uncertainty_var_vec = np.dot(np.transpose(np.square(std_geno_mat[valid_variant_indices, :])), residual_variance_vec[valid_variant_indices])

            # Add to results (one entry per gene-individual pair)
            n_predictive_variants = np.sum(valid_variant_indices)
            for sample_index, sample_name in enumerate(ordered_sample_names):
                gene_individual_key = short_gene_id + ':' + sample_name
                if gene_individual_key in gene_individual_to_results:
                    print('gene-individual repeat assumption erororo')
                    pdb.set_trace()
                gene_individual_to_results[gene_individual_key] = {}
                gene_individual_to_results[gene_individual_key]['observed_expr'] = expr_vec[sample_index]
                gene_individual_to_results[gene_individual_key]['predicted_expr'] = predicted_expr_vec[sample_index]
                gene_individual_to_results[gene_individual_key]['predicted_expr_uncertainty_var'] = predicted_expr_uncertainty_var_vec[sample_index]
                gene_individual_to_results[gene_individual_key]['n_predictive_variants'] = n_predictive_variants

        f.close()

    return gene_individual_to_results


parser = argparse.ArgumentParser(description='Predict expression differences across tissue pairs in shared individuals.')
parser.add_argument('--genotype-stem', dest='genotype_stem', required=True, help='Stem of per-chromosome genotype files.')
parser.add_argument('--sldmc-results-file', dest='sldmc_results_file', required=True, help='SLDMC cross-tissue meta-analyzed results file.')
parser.add_argument('--expr-differences-output-file', dest='expr_differences_output_file', required=True, help='Output file summarizing predicted expression differences.')
parser.add_argument('--anno-method', dest='anno_method', required=True, help='Name of the annotation method used to generate the borzoi annotation files.')
parser.add_argument('--borzoi-effect-file1', dest='borzoi_effect_file1', required=True, help='Borzoi effect sizes file for tissue 1.')
parser.add_argument('--borzoi-annotation-file1', dest='borzoi_annotation_file1', required=True, help='Borzoi variant-gene annotation file for tissue 1.')
parser.add_argument('--genotype-sample-mapping-file1', dest='genotype_sample_mapping_file1', required=True, help='File mapping genotype samples to expression samples for tissue 1.')
parser.add_argument('--expr-file1', dest='expr_file1', required=True, help='Residualized expression file for tissue 1.')
parser.add_argument('--borzoi-effect-file2', dest='borzoi_effect_file2', required=True, help='Borzoi effect sizes file for tissue 2.')
parser.add_argument('--borzoi-annotation-file2', dest='borzoi_annotation_file2', required=True, help='Borzoi variant-gene annotation file for tissue 2.')
parser.add_argument('--genotype-sample-mapping-file2', dest='genotype_sample_mapping_file2', required=True, help='File mapping genotype samples to expression samples for tissue 2.')
parser.add_argument('--expr-file2', dest='expr_file2', required=True, help='Residualized expression file for tissue 2.')
args = parser.parse_args()

genotype_stem = args.genotype_stem
sldmc_results_file = args.sldmc_results_file
expr_differences_output_file = args.expr_differences_output_file
anno_method = args.anno_method
borzoi_effect_file1 = args.borzoi_effect_file1
borzoi_annotation_file1 = args.borzoi_annotation_file1
genotype_sample_mapping_file1 = args.genotype_sample_mapping_file1
expr_file1 = args.expr_file1
borzoi_effect_file2 = args.borzoi_effect_file2
borzoi_annotation_file2 = args.borzoi_annotation_file2
genotype_sample_mapping_file2 = args.genotype_sample_mapping_file2
expr_file2 = args.expr_file2




##################################
# Observed and predicted expression for all gene-individual pairs, in each tissue
##################################
gene_individual_to_results1 = get_observed_and_predicted_expression_for_all_gene_individual_pairs_in_a_tissue(genotype_stem, sldmc_results_file, anno_method, borzoi_effect_file1, borzoi_annotation_file1, genotype_sample_mapping_file1, expr_file1)

gene_individual_to_results2 = get_observed_and_predicted_expression_for_all_gene_individual_pairs_in_a_tissue(genotype_stem, sldmc_results_file, anno_method, borzoi_effect_file2, borzoi_annotation_file2, genotype_sample_mapping_file2, expr_file2)


##################################
# Print results for gene-individual pairs shared between the two tissues
##################################
t = gzip.open(expr_differences_output_file, 'wt')
t.write('gene_id\tindividual_id')
t.write('\tobserved_expr1\tpredicted_expr1\tpredicted_expr_uncertainty_var1')
t.write('\tobserved_expr2\tpredicted_expr2\tpredicted_expr_uncertainty_var2')
t.write('\tobserved_expr_difference\tpredicted_expr_difference\tpredicted_expr_difference_uncertainty_var\n')

n_shared_pairs = 0
for gene_individual_key in [*gene_individual_to_results1]:
    if gene_individual_key not in gene_individual_to_results2:
        continue
    results1 = gene_individual_to_results1[gene_individual_key]
    results2 = gene_individual_to_results2[gene_individual_key]

    observed_expr_difference = results1['observed_expr'] - results2['observed_expr']
    predicted_expr_difference = results1['predicted_expr'] - results2['predicted_expr']
    # Assumes independent prediction errors across the two tissues
    predicted_expr_difference_uncertainty_var = results1['predicted_expr_uncertainty_var'] + results2['predicted_expr_uncertainty_var']

    gene_id, individual_id = gene_individual_key.split(':')
    t.write(gene_id + '\t' + individual_id)
    t.write('\t' + str(results1['observed_expr']) + '\t' + str(results1['predicted_expr']) + '\t' + str(results1['predicted_expr_uncertainty_var']))
    t.write('\t' + str(results2['observed_expr']) + '\t' + str(results2['predicted_expr']) + '\t' + str(results2['predicted_expr_uncertainty_var']))
    t.write('\t' + str(observed_expr_difference) + '\t' + str(predicted_expr_difference) + '\t' + str(predicted_expr_difference_uncertainty_var) + '\n')
    n_shared_pairs = n_shared_pairs + 1

t.close()
print(str(n_shared_pairs) + ' shared gene-individual pairs written to ' + expr_differences_output_file)


