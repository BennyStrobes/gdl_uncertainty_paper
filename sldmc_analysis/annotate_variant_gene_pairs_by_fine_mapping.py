import numpy as np
import sys
import pdb
import gzip
import argparse


# Name of the fine-mapping annotation column (and its two categories)
fine_mapping_anno_name = 'susie_fine_mapped'
fine_mapping_category_names = ['not_fine_mapped', 'fine_mapped']


def normalize_variant_id(variant_id):
	# GTEx hg38 variant ids may or may not carry the '_b38' suffix (chr1_13550_G_A vs chr1_13550_G_A_b38)
	if variant_id.endswith('_b38'):
		return variant_id[:-4]
	return variant_id


def extract_confidently_fine_mapped_vg_pairs(fine_mapping_summary_file, tissue_name, fine_mapping_method, pip_threshold):
	# Returns set of 'variant:gene' names (variant without '_b38' suffix, gene without ensembl version)
	# that are fine-mapped with PIP > pip_threshold in this tissue by this fine-mapping method.
	# Also returns the number of rows found for this (tissue, method), so we can check the tissue exists.
	fine_mapped_vg_pairs = set()
	n_tissue_method_rows = 0

	f = open(fine_mapping_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = np.asarray(data)
			column_indices = {}
			for column_name in ['variant_hg38', 'method', 'tissue', 'gene', 'pip']:
				tmp = np.where(header == column_name)[0]
				if len(tmp) != 1:
					print('assumption eroror: column ' + column_name + ' not found in ' + fine_mapping_summary_file)
					pdb.set_trace()
				column_indices[column_name] = tmp[0]
			continue

		if data[column_indices['method']] != fine_mapping_method:
			continue
		if data[column_indices['tissue']] != tissue_name:
			continue
		n_tissue_method_rows = n_tissue_method_rows + 1

		pip = float(data[column_indices['pip']])
		if pip > 1.0 + 1e-6 or pip < 0.0:
			print('assumption eroror: pip out of range')
			pdb.set_trace()
		if pip <= pip_threshold:
			continue

		variant_id = normalize_variant_id(data[column_indices['variant_hg38']])
		gene_id = data[column_indices['gene']].split('.')[0]
		fine_mapped_vg_pairs.add(variant_id + ':' + gene_id)
	f.close()

	return fine_mapped_vg_pairs, n_tissue_method_rows


def create_annotation_file(borzoi_effect_file, borzoi_annotation_file, fine_mapped_vg_pairs):
	# Same layout as annotate_variant_gene_pairs.py: first 6 columns of the borzoi effect file, then
	# one column per annotation holding the category index of the variant-gene pair.
	# Annotations: intercept (always 0) and the fine-mapping annotation (1 if PIP > threshold, else 0).
	category_counts = np.zeros(len(fine_mapping_category_names), dtype=int)
	matched_fine_mapped_vg_pairs = set()

	f = gzip.open(borzoi_effect_file,'rt')
	t = gzip.open(borzoi_annotation_file, 'wt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
			t.write('\tintercept\t' + fine_mapping_anno_name + '\n')
			continue

		if data[4] == data[5]:
			print('assumption eroroor')
			pdb.set_trace()

		# Borzoi effect file: gene (no ensembl version), variant, chr, snp_pos, a0, a1, borzoi_effect_size
		vg_name = normalize_variant_id(data[1]) + ':' + data[0].split('.')[0]
		if vg_name in fine_mapped_vg_pairs:
			category_index = 1
			matched_fine_mapped_vg_pairs.add(vg_name)
		else:
			category_index = 0
		category_counts[category_index] = category_counts[category_index] + 1

		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5])
		t.write('\t0\t' + str(category_index) + '\n')

	f.close()
	t.close()

	return category_counts, matched_fine_mapped_vg_pairs


def create_annotation_category_file(annotation_category_file):
	t = open(annotation_category_file, 'w')
	t.write('anno_name\tsource\tcategory_index\tcategory_name\n')
	t.write('intercept\tcustom\t0\tintercept\n')
	for category_iter, category_name in enumerate(fine_mapping_category_names):
		t.write(fine_mapping_anno_name + '\tfine_mapping\t' + str(category_iter) + '\t' + category_name + '\n')
	t.close()
	return


def create_annotation_category_counts_file(annotation_category_counts_file, category_counts):
	t = open(annotation_category_counts_file, 'w')
	t.write('anno_name\tsource\tcategory_index\tcategory_name\tn_variant_gene_pairs\n')
	t.write('intercept\tcustom\t0\tintercept\t' + str(np.sum(category_counts)) + '\n')
	t.write('intercept\tcustom\t-1\tno_category\t0\n')
	for category_iter, category_name in enumerate(fine_mapping_category_names):
		t.write(fine_mapping_anno_name + '\tfine_mapping\t' + str(category_iter) + '\t' + category_name + '\t' + str(category_counts[category_iter]) + '\n')
	t.write(fine_mapping_anno_name + '\tfine_mapping\t-1\tno_category\t0\n')
	t.close()
	return


#####################
# Command line args
#####################
parser = argparse.ArgumentParser(description='Annotate variant-gene pairs by whether they are confidently fine-mapped (PIP > threshold) in a given tissue.')
parser.add_argument('--borzoi-effect-file', dest='borzoi_effect_file', required=True, help='Borzoi effect sizes file (the variant-gene pairs to annotate).')
parser.add_argument('--fine-mapping-summary-file', dest='fine_mapping_summary_file', required=True, help='GTEx fine-mapping summary file (GTEx_49tissues_release1.tsv).')
parser.add_argument('--tissue-name', dest='tissue_name', required=True, help='GTEx tissue name, as it appears in the tissue column of the fine-mapping summary file.')
parser.add_argument('--borzoi-annotation-file', dest='borzoi_annotation_file', required=True, help='Output variant-gene annotation file. The (annotation, category) files are written alongside it.')
parser.add_argument('--fine-mapping-method', dest='fine_mapping_method', default='SUSIE', help='Fine-mapping method to use (method column of the fine-mapping summary file).')
parser.add_argument('--pip-threshold', dest='pip_threshold', type=float, default=0.9, help='Variant-gene pairs with PIP strictly greater than this are fine-mapped.')
args = parser.parse_args()


# Describes the (annotation, category) pairs making up the annotation file
annotation_category_file = args.borzoi_annotation_file.split('.txt.gz')[0] + '_categories.txt'
# Same rows plus the number of variant-gene pairs falling in each category
annotation_category_counts_file = args.borzoi_annotation_file.split('.txt.gz')[0] + '_category_counts.txt'


# Extract confidently fine-mapped variant-gene pairs in this tissue
fine_mapped_vg_pairs, n_tissue_method_rows = extract_confidently_fine_mapped_vg_pairs(args.fine_mapping_summary_file, args.tissue_name, args.fine_mapping_method, args.pip_threshold)
if n_tissue_method_rows == 0:
	print('assumption eroror: no ' + args.fine_mapping_method + ' rows for tissue ' + args.tissue_name + ' in ' + args.fine_mapping_summary_file)
	pdb.set_trace()
print(str(n_tissue_method_rows) + ' ' + args.fine_mapping_method + ' variant-gene pairs in ' + args.tissue_name)
print(str(len(fine_mapped_vg_pairs)) + ' of them with PIP > ' + str(args.pip_threshold))

# Print annotation file
category_counts, matched_fine_mapped_vg_pairs = create_annotation_file(args.borzoi_effect_file, args.borzoi_annotation_file, fine_mapped_vg_pairs)
print(str(np.sum(category_counts)) + ' variant-gene pairs annotated')
print('\t' + fine_mapping_category_names[0] + ': ' + str(category_counts[0]))
print('\t' + fine_mapping_category_names[1] + ': ' + str(category_counts[1]))
print(str(len(matched_fine_mapped_vg_pairs)) + ' of ' + str(len(fine_mapped_vg_pairs)) + ' fine-mapped variant-gene pairs found in borzoi effect file')
if category_counts[1] == 0:
	print('WARNING: no variant-gene pairs fall in the fine_mapped category')

# Print (annotation, category) pair file
create_annotation_category_file(annotation_category_file)

# Print (annotation, category) pair counts file
create_annotation_category_counts_file(annotation_category_counts_file, category_counts)
