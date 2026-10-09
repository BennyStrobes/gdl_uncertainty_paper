import numpy as np
import sys
import pdb
import gzip


def create_mapping_from_vg_pair_to_pip(susie_fine_mapping_file):
	f = gzip.open(susie_fine_mapping_file, 'rt')
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene = data[0]
		variant = data[1]
		pip = float(data[2])
		vg_pair = variant + ':' + gene
		if vg_pair in mapping:
			print('repeat variant-gene pair assumption error')
			pdb.set_trace()
		mapping[vg_pair] = pip
	f.close()
	return mapping


##########################
# Command line args
##########################
sldmc_variant_gene_annotation_file = sys.argv[1]  # Existing SLDMC annotation file (defines universe of variant-gene pairs)
susie_fine_mapping_file = sys.argv[2]
pip_thresh = float(sys.argv[3])
fm_status_sldmc_annotation_file = sys.argv[4]
fm_status_sldmc_annotation_category_file = sys.argv[5]  # Must be annotation file name with '.txt.gz' swapped for '_categories.txt'


# Create mapping from variant-gene pair to PIP
vg_to_pip = create_mapping_from_vg_pair_to_pip(susie_fine_mapping_file)


# Write companion category file
# Two categories: 0 = not confidently fine-mapped (PIP < thresh, or missing from susie file), 1 = confidently fine-mapped (PIP >= thresh)
t_cat = open(fm_status_sldmc_annotation_category_file, 'w')
t_cat.write('anno_name\tsource\tcategory_index\tcategory_name\n')
t_cat.write('fm_status\tsusie_pip_' + str(pip_thresh) + '\t0\tnot_fine_mapped\n')
t_cat.write('fm_status\tsusie_pip_' + str(pip_thresh) + '\t1\tfine_mapped\n')
t_cat.close()


# Write annotation file: six shared columns from existing SLDMC annotation file, then fm_status
f = gzip.open(sldmc_variant_gene_annotation_file, 'rt')
t = gzip.open(fm_status_sldmc_annotation_file, 'wt')
t.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\tfm_status\n')
head_count = 0
n_pairs = 0
n_fine_mapped = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene = data[0]
	variant = data[1]
	vg_pair = variant + ':' + gene
	fm_status = 0
	if vg_pair in vg_to_pip and vg_to_pip[vg_pair] >= pip_thresh:
		fm_status = 1
	n_pairs = n_pairs + 1
	n_fine_mapped = n_fine_mapped + fm_status
	t.write('\t'.join(data[:6]) + '\t' + str(fm_status) + '\n')
f.close()
t.close()

print(str(n_fine_mapped) + ' of ' + str(n_pairs) + ' variant-gene pairs confidently fine-mapped (PIP >= ' + str(pip_thresh) + ')')
print(fm_status_sldmc_annotation_file)
