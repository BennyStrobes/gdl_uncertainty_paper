import numpy as np
import os
import sys
import pdb
import h5py
import gzip
import glob
import pyarrow.parquet as pq



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
	file_name = os.path.basename(borzoi_input_file)
	return int(file_name.split('_chunk_')[1].split('_')[0])


def stream_borzoi_results(borzoi_results_dir, borzoi_target_index, chunk_size=50000):
	file_pattern = os.path.join(borzoi_results_dir, 'model_0_chunk_*_borzoi_gtex_only_results.h5')
	borzoi_input_files = sorted(glob.glob(file_pattern), key=extract_chunk_num_from_borzoi_file)

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








def load_in_eqtl_vg_pairs(filer):
	dicti = {}
	f = gzip.open(filer,'rt')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[1]
		gene_id = data[0]
		if variant_id + ':' + gene_id in dicti:
			print('erroro')
			pdb.set_trace()
		dicti[variant_id + ':' + gene_id] = 1
	f.close()


	return dicti


def extract_info_on_each_variant_gene_pair(eqtl_sumstats_dir, tissue_name):
	dicti = {}
	for chrom_num in range(1,23):
		filer = eqtl_sumstats_dir + tissue_name + '.v10.allpairs.chr' + str(chrom_num) + '.parquet'
		pf = pq.ParquetFile(filer)

		for rg in range(pf.num_row_groups):
			table = pf.read_row_group(
				rg,
				columns=['gene_id', 'variant_id', 'tss_distance', 'af', 'slope', 'slope_se', 'ma_count']
			)

			if table.num_columns == 0:
				continue

			gene_ids = table['gene_id'].to_pylist()
			variant_ids = table['variant_id'].to_pylist()
			tss_dists = table['tss_distance'].to_pylist()
			afs = table['af'].to_pylist()

			for ii, gene_id in enumerate(gene_ids):
				maf = afs[ii]
				if maf > 0.5:
					maf = 1.0 - maf
				if maf == 0.0:
					continue

				abs_tss_dist = np.abs(tss_dists[ii])
				if abs_tss_dist > 100000:
					continue

				variant_id = variant_ids[ii]
				vg_name = variant_id + ':' + gene_id.split('.')[0]
				dicti[vg_name] = maf
	return dicti

###################
# Command line args
####################
gtex_v10_pc_genes_gtf = sys.argv[1]
borzoi_results_dir = sys.argv[2]
borzoi_target_index = int(sys.argv[3])
target_tissue = sys.argv[4]
target_sample = sys.argv[5]
borzoi_output_dir = sys.argv[6]


# Extract dictionary list of protein coding genes
pc_genes = extract_dictionary_list_of_protein_coding_genes(gtex_v10_pc_genes_gtf)


# Open output file handle
output_file = borzoi_output_dir + target_tissue + '_' + target_sample + '_borzoi_effects.txt.gz'
t = gzip.open(output_file,'wt')
t.write('gene\tvariant\tchr\tsnp_pos\ta0\ta1\tborzoi_effect_size\n')


for snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects in stream_borzoi_results(borzoi_results_dir, borzoi_target_index):
	for chrom, pos, variant_id, gene_id, borzoi_effect_size in zip(snp_chrom, snp_pos, snp_ids, gene_ids, borzoi_effects):
		gene_name = gene_id.split('.')[0]
		if gene_name not in pc_genes:
			continue

		var_info = variant_id.split('_')
		if len(var_info) < 4:
			continue
		if len(var_info[2]) != 1 or len(var_info[3]) != 1:
			continue

		chrom_string = str(chrom)
		if chrom_string.startswith('chr'):
			chrom_string = chrom_string.split('hr')[1]


		t.write(gene_name + '\t' + variant_id + '\t' + chrom_string + '\t' + str(pos) + '\t' + var_info[2] + '\t' + var_info[3] + '\t' + str(borzoi_effect_size) + '\n')

t.close()
