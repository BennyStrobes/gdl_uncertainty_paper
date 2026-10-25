import numpy as np
import os
import sys
import pdb
import h5py
import gzip
import glob



def extract_dictionary_list_of_protein_coding_genes(pc_genes_gtf):
	valid_chroms = {}
	for chrom_num in range(1,23):
		valid_chroms['chr' + str(chrom_num)] = 1


	f = open(pc_genes_gtf)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if data[0] not in valid_chroms:
			continue
		ens_id = data[8].split(';')[0].split('"')[1]
		if ens_id.startswith('ENSG') == False:
			print('assumption oernroro')
			pdb.set_trace()
		dicti[ens_id.split('.')[0]] = 1

	f.close()

	return dicti


def decode_if_bytes(arr):
	if len(arr) == 0:
		return arr
	if isinstance(arr[0], bytes):
		return np.asarray([x.decode('utf-8') for x in arr])
	return arr


def extract_chunk_num_from_borzoi_file(borzoi_input_file):
	# File names like: variants_chunk_04.borzoi_output.h5
	file_name = os.path.basename(borzoi_input_file)
	return int(file_name.split('variants_chunk_')[1].split('.')[0])


def get_ordered_borzoi_input_files(borzoi_results_dir):
	file_pattern = os.path.join(borzoi_results_dir, 'variants_chunk_*.borzoi_output.h5')
	borzoi_input_files = sorted(glob.glob(file_pattern), key=extract_chunk_num_from_borzoi_file)
	if len(borzoi_input_files) == 0:
		print('assumption eroror: no borzoi h5 files found matching ' + file_pattern)
		pdb.set_trace()
	return borzoi_input_files


def count_lines_minus_header(filer):
	f = open(filer)
	n_lines = 0
	for line in f:
		n_lines = n_lines + 1
	f.close()
	return n_lines - 1


def select_borzoi_target_column(borzoi_input_files, orig_target_index, gtex_only_target_index, n_gtex_only_targets):
	# The h5 files have the same format as the genome-wide borzoi run (one logRef/logAlt column per target),
	# but may have been run on either the full borzoi target set (index with orig_target_index) or on the
	# gtex-only target set (index with gtex_only_target_index). Pick based on the number of target columns.
	with h5py.File(borzoi_input_files[0], 'r') as f:
		n_targets = f['logRef'].shape[1]
	if n_targets == n_gtex_only_targets:
		print('h5 files have ' + str(n_targets) + ' targets (gtex-only target set): using gtex-only target index ' + str(gtex_only_target_index))
		return gtex_only_target_index
	if orig_target_index >= n_targets:
		print('assumption eroror: h5 files have ' + str(n_targets) + ' targets, which matches neither the gtex-only target set (' + str(n_gtex_only_targets) + ') nor the full target set')
		pdb.set_trace()
	print('h5 files have ' + str(n_targets) + ' targets (assumed full target set): using original target index ' + str(orig_target_index))
	return orig_target_index


def stream_borzoi_results(borzoi_input_files, borzoi_target_index, chunk_size=50000):
	for borzoi_input_file in borzoi_input_files:
		with h5py.File(borzoi_input_file, 'r') as f:
			snp_chrom_ds = f['snp_chrom']
			snp_pos_ds = f['snp_pos']
			snp_ds = f['snp']
			gene_ds = f['gene']
			logRef_ds = f['logRef']
			logAlt_ds = f['logAlt']

			n_rows = gene_ds.shape[0]
			for start in range(0, n_rows, chunk_size):
				end = min(start + chunk_size, n_rows)

				snp_chrom = decode_if_bytes(snp_chrom_ds[start:end])
				snp_pos = snp_pos_ds[start:end]
				snp_ids = decode_if_bytes(snp_ds[start:end])
				gene_ids = decode_if_bytes(gene_ds[start:end])
				logRef = logRef_ds[start:end, borzoi_target_index]
				logAlt = logAlt_ds[start:end, borzoi_target_index]
				borzoi_effects = logAlt - logRef

				yield snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects



###################
# Command line args
####################
gtex_v10_pc_genes_gtf = sys.argv[1]
borzoi_results_dir = sys.argv[2]
orig_target_index = int(sys.argv[3])
gtex_only_target_index = int(sys.argv[4])
borzoi_gtex_only_target_names_file = sys.argv[5]
target_tissue = sys.argv[6]
target_sample = sys.argv[7]
borzoi_output_dir = sys.argv[8]


# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)

# Borzoi h5 files and the logRef/logAlt column corresponding to this target
borzoi_input_files = get_ordered_borzoi_input_files(borzoi_results_dir)
n_gtex_only_targets = count_lines_minus_header(borzoi_gtex_only_target_names_file)
borzoi_target_index = select_borzoi_target_column(borzoi_input_files, orig_target_index, gtex_only_target_index, n_gtex_only_targets)


# Open output file handle
output_file = borzoi_output_dir + target_tissue + '_' + target_sample + '_rare_variant_borzoi_effects.txt.gz'
t = gzip.open(output_file,'wt')
t.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\tborzoi_effect_size\n')

n_written = 0
n_unparsed_ids = 0
for snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects in stream_borzoi_results(borzoi_input_files, borzoi_target_index):
	for chrom, pos, variant_id, gene_id, borzoi_effect_size in zip(snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects):
		gene_name = gene_id.split('.')[0]
		if gene_name not in pc_genes:
			continue

		# Variant ids like chr1_12345_A_G_b38 (alleles are taken from the id)
		var_info = variant_id.split('_')
		if len(var_info) < 4:
			n_unparsed_ids = n_unparsed_ids + 1
			continue
		if len(var_info[2]) != 1 or len(var_info[3]) != 1:
			continue
		if var_info[1] != str(pos):
			print('assumption eroror: variant id position does not match snp_pos: ' + variant_id + ' ' + str(pos))
			pdb.set_trace()

		chrom_string = str(chrom)
		if chrom_string.startswith('chr'):
			chrom_string = chrom_string.split('hr')[1]

		# Rewrite the variant id in the format of the rare variant carrier tables (chr1_12345_A_G_b38)
		standardized_variant_id = 'chr' + chrom_string + '_' + str(pos) + '_' + var_info[2] + '_' + var_info[3] + '_b38'

		t.write(gene_name + '\t' + standardized_variant_id + '\t' + chrom_string + '\t' + str(pos) + '\t' + var_info[2] + '\t' + var_info[3] + '\t' + str(borzoi_effect_size) + '\n')
		n_written = n_written + 1

t.close()

print(str(n_written) + ' variant-gene pairs written')
if n_unparsed_ids > 0:
	print('WARNING: ' + str(n_unparsed_ids) + ' variant-gene pairs skipped because the variant id could not be parsed into chrom_pos_ref_alt')
if n_written == 0:
	print('assumption eroror: no variant-gene pairs written')
	pdb.set_trace()
