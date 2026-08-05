import numpy as np
import sys
import pdb
import gzip


def create_mapping_from_variant_id_to_genotype_std(est_eqtl_effect_size_file):
    """
    Create mapping from variant id to genotype standard deviation
    :param est_eqtl_effect_size_file: file containing estimated eqtl effect sizes
    :return: mapping from variant id to genotype standard deviation
    """
    variant_id_to_genotype_std = {}
    with gzip.open(est_eqtl_effect_size_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            variant_id = fields[header.index('variant')]
            genotype_std = float(fields[header.index('genotype_sdev')])
            if variant_id in variant_id_to_genotype_std:
                if variant_id_to_genotype_std[variant_id] != genotype_std:
                    print(f"Warning: variant {variant_id} has different genotype std values: {variant_id_to_genotype_std[variant_id]} and {genotype_std}") 
                    pdb.set_trace()
            variant_id_to_genotype_std[variant_id] = genotype_std
    return variant_id_to_genotype_std









######################
# Command line args
#######################
est_eqtl_effect_size_file = sys.argv[1]
est_borzoi_standardized_effect_size_file = sys.argv[2]
est_borzoi_effect_size_file = sys.argv[3]

# First create mappping from variant id to genotype standard deviation
variant_id_to_genotype_std = create_mapping_from_variant_id_to_genotype_std(est_eqtl_effect_size_file)

# Now scale the borzoi standardized effect sizes to per-allele effect sizes
# and write to output file
with gzip.open(est_borzoi_standardized_effect_size_file, 'rt') as f_in, gzip.open(est_borzoi_effect_size_file, 'wt') as f_out:
    header = f_in.readline().strip().split('\t')
    f_out.write('\t'.join(header) + '\n')
    for line in f_in:
        fields = line.strip().split('\t')
        variant_id = fields[header.index('variant')]
        borzoi_standardized_effect = float(fields[header.index('borzoi_effect_size')])
        if variant_id not in variant_id_to_genotype_std:
            print(f"Warning: variant {variant_id} not found in genotype std mapping.")
            borzoi_per_allele_effect = 0
        else:
            genotype_std = variant_id_to_genotype_std[variant_id]
            if genotype_std == 0:
                borzoi_per_allele_effect = 0 # Ok. because gets multiplied by genotype std in downstream analysis
            else:
               borzoi_per_allele_effect = borzoi_standardized_effect / genotype_std
        fields[header.index('borzoi_effect_size')] = str(borzoi_per_allele_effect)
        f_out.write('\t'.join(fields) + '\n')