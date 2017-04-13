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

def parse_pop_map(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, sara_id, origin, country in izip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country}
    return pop_map

def kinship_all_genes(snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5',
                 min_maf=0.05,
                 max_strain_num=198):
    """
    Calculates the kinship. Our expectation is that will be close to identity matrix
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

    tuple_index_kinship = (K_snps, strain_mask)

    headers = list()
    maps = parse_pop_map()
    for i in ordered_strains:
        headers.append( 'SM' +maps[i]['sara_id'])

    np.savetxt('strains_order.csv', headers,  delimiter=",", fmt='%s')
    np.savetxt('kinship.csv', K_snps, fmt='%.18e', delimiter=',', header = str(headers))
    # Saving as a compressed numpy file:
    #np.savez_compressed("{}/{}".format(kinship_matrix), matrix=snps, strains=strains, maf = maf) 

    return(tuple_index_kinship)
#kinship_all_genes()


def kinship_pseudo_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/',
                        num_strains = 198):
    # Changing directory
    os.chdir(directory)

    genes = []
    for f in glob.glob('*.npz'):
        with np.load(directory + f) as data:
            genes.append((f, data["matrix"], data["strains"]))

    K_snps = sp.zeros((num_strains, num_strains))
    counts_mat_snps = sp.zeros((num_strains, num_strains))

    for gene in genes:
        name, snps, strain_mask = (gene)
        K_snps_slice = K_snps[strain_mask]
        K_snps_slice[:, strain_mask] += sp.dot(snps, snps.T)
        K_snps[strain_mask] = K_snps_slice
        counts_mat_snps_slice = counts_mat_snps[strain_mask]
        counts_mat_snps_slice[:, strain_mask] += len(snps)
        counts_mat_snps[strain_mask] = counts_mat_snps_slice

        if len(strain_mask) == 198:
            final_strain_mask = strain_mask 

    K_snps = K_snps / counts_mat_snps  # element-wise division
    print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))

    pl.matshow(K_snps)
    pl.title('Kinship pseudo genes')
    pl.colorbar()
    pl.savefig('kinshi_pseudo_genes.png')
    #pl.show()

    tuple_index_kinship = (K_snps, final_strain_mask)

    return(tuple_index_kinship)

#print kinship_pseudo_genes()

def kinship_versus_corrected_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/'):
    os.chdir(directory)

	# Upload overall kinship matrix
    kinship, k_strains = (kinship_pseudo_genes())

	#This is the gene SNPs matrix
    genes = []
    for f in glob.glob('*.npz'):
        with np.load(directory + f) as data:
	    	genes.append((f, data["matrix"], data["strains"]))

    r_scores = []
    gene_name = []
    print len(genes)
    for gene in genes:

        name, snps, strains_1 = (gene)
        strains_mask_1 = np.in1d(strains_1, k_strains, assume_unique = True)
        filtered_strains_1 = strains_1[strains_mask_1]

        strain_mask_2 = np.in1d(k_strains, filtered_strains_1, assume_unique=True)
        filtered_strains_2 = k_strains[strain_mask_2]

		# Construct GRM for a singular gene
        total_snps_1 = snps[strains_mask_1, :]
        grm_1 = np.divide(np.dot(snps, snps.T), snps.shape[1])

        flat_grm_1 = grm_1.flatten()
        norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
        norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))
		
		# Overall kinship
        grm_2 = kinship[strain_mask_2, :]
        grm_2 = grm_2[:, strain_mask_2]

        flat_grm_2 = grm_2.flatten()
        norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
        norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))

        gene_name.append(gene[0][:-4])
        r_scores.append(pearsonr(norm_flat_grm1, norm_flat_grm2)[0])

    # Include number of markers for each gene
    # Include the plasmid origin
    # Include the gene functionality
    LD_stats = pd.DataFrame({'r_scores': r_scores,'gene':gene_name})
    LD_stats.to_csv('introgressed_gene_stats.csv', header = True)

kinship_versus_corrected_genes()