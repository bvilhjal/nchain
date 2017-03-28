# Identify genes that do not follow the overall population structure, i. e. genes that introgress.
# For this purpose we need to:
# 1. Construct the phylogenetic tree based on all the genes
# 2. Use the correct population structure genes and correlates it to the kinship matrix
import h5py
import numpy as np
import scipy as sp
import pandas as pd
import os
import glob 
from scipy.stats.stats import pearsonr
import pylab as pl

def kinship_all_genes(snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5',
                 min_maf=0.1,
                 max_strain_num=198):
    """
    Calculates the kinship
    """
    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    all_strains = set()
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) == max_strain_num:
            all_strains = set(strains).union(all_strains)
    num_strains = len(all_strains)
    print 'Found %d "distinct" strains' % num_strains

    ordered_strains = sorted(list(all_strains))

    strain_index = pd.Index(ordered_strains)
    K_snps = sp.zeros((num_strains, num_strains))
    counts_mat_snps = sp.zeros((num_strains, num_strains))

    for i, gg in enumerate(gene_groups):
        if i % 100 == 0:
            print 'Working on gene nr. %d' % i
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) < max_strain_num:
            strain_mask = strain_index.get_indexer(strains)

            # Already normalized snps
            snps = data_g['norm_snps'][...]
            freqs = data_g['freqs'][...]
            mafs = sp.minimum(freqs, 1 - freqs)
            maf_mask = mafs > min_maf
            snps = snps[maf_mask]

            if len(snps) == 0:
                continue
            K_snps_slice = K_snps[strain_mask]
            K_snps_slice[:, strain_mask] += sp.dot(snps.T, snps)
            K_snps[strain_mask] = K_snps_slice
            counts_mat_snps_slice = counts_mat_snps[strain_mask]
            counts_mat_snps_slice[:, strain_mask] += len(snps)
            counts_mat_snps[strain_mask] = counts_mat_snps_slice

    K_snps = K_snps / counts_mat_snps  # element-wise division
    print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))
    #correlation_plot(K_snps, wrt=False)

    tuple_index_kinship = (K_snps, strain_mask)

    #np.savez_compressed("{}/{}".format(kinship_matrix), matrix=snps, strains=strains, maf = maf) 
    return(tuple_index_kinship)
#kinship_all_genes()

def kinship_versus_genes_wo_correction(snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5',
										max_strain_num = 198):
    
	# Upload overall kinship matrix
	kinship, k_strains = (kinship_all_genes())
	
	h5f = h5py.File(snps_file)
	gene_groups = h5f.keys()
	all_strains = set()
	for gg in gene_groups:
		data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) == max_strain_num:
            all_strains = set(strains).union(all_strains)
    
	num_strains = len(all_strains)
	print 'Found %d "distinct" strains' % num_strains

	ordered_strains = sorted(list(all_strains))

	strain_index = pd.Index(ordered_strains)
	K_snps = sp.zeros((num_strains, num_strains))
	counts_mat_snps = sp.zeros((num_strains, num_strains))

	for i, gg in enumerate(gene_groups):
		if i % 100 == 0:
			print 'Working on gene nr. %d' % i
    	data_g = h5f[gg]
    	strains = data_g['strains'][...]
    	if len(strains) < max_strain_num:
        	strain_mask = strain_index.get_indexer(strains)

            # Already normalized snps
        	snps = data_g['norm_snps'][...]
        	freqs = data_g['freqs'][...]
        	mafs = sp.minimum(freqs, 1 - freqs)
        	maf_mask = mafs > min_maf
        	snps = snps[maf_mask]

#kinship_versus_genes_wo_correction()


def kinship_versus_corrected_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/'):
	os.chdir(directory)

	# Upload overall kinship matrix
	kinship, k_strains = (kinship_all_genes())

	#This is the gene SNPs matrix
	genes = []
	for f in glob.glob('*.npz'):
	    with np.load(directory + f) as data:
	    	genes.append((f, data["matrix"], data["strains"], data["maf"]))

	r_scores = []
	gene_name = []
	print len(genes)
	for gene in genes:

		name, snps, strains_1, maf = (gene)

		strains_mask_1 = np.in1d(strains_1, k_strains, assume_unique = True)
		filtered_strains_1 = strains_1[strains_mask_1]

		strain_mask_2 = np.in1d(k_strains, filtered_strains_1, assume_unique=True)
		filtered_strains_2 = k_strains[strain_mask_2]

		# Construct GRM for a singular gene
		total_snps_1 = snps[strains_mask_1, :]
		grm_1 = np.divide(np.dot(snps, snps.T), snps.shape[1])
		pl.matshow(grm_1)
		pl.show()

		flat_grm_1 = grm_1.flatten()
		norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
		norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))
		
		# Overall kinship
		grm_2 = kinship[strain_mask_2, :]
		grm_2 = grm_2[:, strain_mask_2]
		pl.matshow(grm_2)
		pl.show()
		flat_grm_2 = grm_2.flatten()
		norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
		norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))

		gene_name.append(gene)
		r_scores.append(pearsonr(norm_flat_grm1, norm_flat_grm2))

	LD_stats = pd.DataFrame({'r_scores': r_scores,
							'gene':gene_name})
	LD_stats.to_scv('introgressed_gene_stats.csv', heder = True)

kinship_versus_corrected_genes()

