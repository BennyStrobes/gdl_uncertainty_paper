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


def extract_annotation_categories(annotation_category_file):
	# Companion file to an SLDMC annotation file. Columns: anno_name  source  category_index  category_name
	# Returns, per annotation, the ordered list of its category names
	anno_name_to_category_names = {}
	f = open(annotation_category_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_name = data[0]
		category_index = int(data[2])
		category_name = data[3]
		if anno_name not in anno_name_to_category_names:
			anno_name_to_category_names[anno_name] = []
		if category_index != len(anno_name_to_category_names[anno_name]):
			print('assumption error: categories not in index order for ' + anno_name)
			pdb.set_trace()
		anno_name_to_category_names[anno_name].append(category_name)
	f.close()
	return anno_name_to_category_names


##########################
# Command line args
##########################
sldmc_variant_gene_annotation_file = sys.argv[1]  # Existing SLDMC annotation file (defines universe of variant-gene pairs and simulated annotations)
susie_fine_mapping_file = sys.argv[2]
pip_thresh = float(sys.argv[3])
fm_status_sldmc_annotation_file = sys.argv[4]
fm_status_sldmc_annotation_category_file = sys.argv[5]  # Must be annotation file name with '.txt.gz' swapped for '_categories.txt'

# Categories of the existing SLDMC annotation file (SLDMC naming convention)
sldmc_annotation_category_file = sldmc_variant_gene_annotation_file.split('.txt.gz')[0] + '_categories.txt'
input_anno_name_to_category_names = extract_annotation_categories(sldmc_annotation_category_file)

# Create mapping from variant-gene pair to PIP
vg_to_pip = create_mapping_from_vg_pair_to_pip(susie_fine_mapping_file)

fm_status_names = ['not_fine_mapped', 'fine_mapped']
source_name = 'susie_pip_' + str(pip_thresh)


##########################
# Pass through existing SLDMC annotation file
##########################
f = gzip.open(sldmc_variant_gene_annotation_file, 'rt')
t = gzip.open(fm_status_sldmc_annotation_file, 'wt')
head_count = 0
n_pairs = 0
n_fine_mapped = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		input_anno_names = data[6:]
		# Annotations to cross with fm_status: every input annotation except the intercept
		crossed_anno_indices = [ii for ii, anno_name in enumerate(input_anno_names) if anno_name != 'intercept']
		for anno_index in crossed_anno_indices:
			if input_anno_names[anno_index] not in input_anno_name_to_category_names:
				print('assumption error: ' + input_anno_names[anno_index] + ' missing from category file')
				pdb.set_trace()
		# Header
		# Output annotations: fm_status (binary), then each input annotation crossed with fm_status
		output_anno_names = ['fm_status'] + [input_anno_names[anno_index] + '_x_fm_status' for anno_index in crossed_anno_indices]
		t.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\t' + '\t'.join(output_anno_names) + '\n')
		continue
	gene = data[0]
	variant = data[1]
	vg_pair = variant + ':' + gene
	fm_status = 0
	if vg_pair in vg_to_pip and vg_to_pip[vg_pair] >= pip_thresh:
		fm_status = 1
	n_pairs = n_pairs + 1
	n_fine_mapped = n_fine_mapped + fm_status

	output_annos = [fm_status]
	for anno_index in crossed_anno_indices:
		input_category_index = int(data[6 + anno_index])
		if input_category_index < 0:
			# Pair in no category of this annotation -> no category of the crossed annotation
			output_annos.append(-1)
		else:
			# Crossed category index: (input category, fm_status) pairs, fm_status fastest-varying
			output_annos.append(input_category_index*2 + fm_status)
	t.write('\t'.join(data[:6]) + '\t' + '\t'.join([str(x) for x in output_annos]) + '\n')
f.close()
t.close()


##########################
# Write companion category file
##########################
t_cat = open(fm_status_sldmc_annotation_category_file, 'w')
t_cat.write('anno_name\tsource\tcategory_index\tcategory_name\n')
for fm_status, fm_status_name in enumerate(fm_status_names):
	t_cat.write('fm_status\t' + source_name + '\t' + str(fm_status) + '\t' + fm_status_name + '\n')
for anno_index in crossed_anno_indices:
	input_anno_name = input_anno_names[anno_index]
	for input_category_index, input_category_name in enumerate(input_anno_name_to_category_names[input_anno_name]):
		for fm_status, fm_status_name in enumerate(fm_status_names):
			t_cat.write(input_anno_name + '_x_fm_status\t' + source_name + '\t' + str(input_category_index*2 + fm_status) + '\t' + input_category_name + '_' + fm_status_name + '\n')
t_cat.close()

print(str(n_fine_mapped) + ' of ' + str(n_pairs) + ' variant-gene pairs confidently fine-mapped (PIP >= ' + str(pip_thresh) + ')')
print(fm_status_sldmc_annotation_file)
