'''Implement a simpler version of the Mantel test, ignoring population structure.'''

# In this script we calculate the kinship matrix for each core gene and do a simple mantel test (correlation of matrix) between pairs of genes

import numpy as np
import scipy as sp
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import seaborn as sns
import time
from collections import OrderedDict

def parse_nod():
    nod_genes = OrderedDict()
    # Full matrix:
    #nod_genes = {4144: 'nodX', 4143: 'nodN', 4142: 'nodM', 4141: 'nodL',10151: 'nodR', 4140: 'nodE', 4139: 'nodF', 4138: 'nodD', 4137: 'nodA', 8237: 'nodB', 4136: 'nodC', 4135: 'nodI', 4134: 'nodJ', 4131: 'fixT', 4130:'fixN', 4129: 'nifB', 4128: 'nifA', 4127: 'fixX', 4126:'fixC', 4125: 'fixB', 4124: 'fixA', 4123: 'nifH', 4122: 'nifD', 4121: 'nifK', 4120: 'nifE', 4119: 'nifN', 2448: 'rpoB', 2140: 'recA'}
    nod_genes = OrderedDict([(4144, 'nodX'), (4143, 'nodN'), (4142, 'nodM'), (4141, 'nodL'), (4140, 'nodE'), (4139, 'nodF'), (4138, 'nodD'), (4137, 'nodA'), (4136, 'nodC'), (4135, 'nodI'), (4134, 'nodJ'), (4129, 'nifB'), (4128, 'nifA'), (4127, 'fixX'), (4126,'fixC'), (4125, 'fixB'), (4124, 'fixA'), (4123, 'nifH'), (4122, 'nifD'), (4121, 'nifK'), (4120, 'nifE'), (2448, 'rpoB'), (2140, 'recA')])
    print nod_genes
    return(nod_genes)

def minor_allele_filter(gene_matrix, maf):
    """Filtering for minor allele frequency, it assumes that the matrix is binary and that is n x m, where columns are markers (m)."""

    freqs = sp.mean(gene_matrix, 0)
    mafs = sp.minimum(freqs, 1 - freqs)

    maf_filter = mafs > maf
    matrix_mafs = gene_matrix[:, maf_filter]

    # Normalizing the columns:
    norm_matrix = (matrix_mafs - np.mean(matrix_mafs, axis=0)) / np.std(matrix_mafs, axis=0)
    return(norm_matrix)

def average_genotype_matrix(X):
    'Computes the average genotype matrix, it assumes that the input matrix (Markers/M in rows and individuals/N in collumns)'
    average_X = (X - np.mean(X, axis = 0))
    return(average_X)

def correlation_plot(df, wrt = True):
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(12, 9))

    # Clustering with seaborn
    with sns.axes_style("white"):
        ax = sns.heatmap(df, square=True, annot= wrt, annot_kws={"size": 9}, cmap="RdYlGn", vmin=0, vmax=1)
    f.tight_layout()
    plt.show()

def kinship_all_genes(snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5',
                 min_maf = 0.05,
                 max_strain_num=200):
    """
    Calculates the kinship
    """
    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    all_strains = set()
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains)<max_strain_num:
            all_strains = set(strains).union(all_strains)
    num_strains = len(all_strains)
    print 'Found %d "distinct" strains'%num_strains
    
    ordered_strains = sorted(list(all_strains))
    
    strain_index = pd.Index(ordered_strains)
    K_snps = sp.zeros((num_strains,num_strains))
    counts_mat_snps = sp.zeros((num_strains,num_strains))
   
    for i, gg in enumerate(gene_groups):
        if i%100==0:
            print 'Working on gene nr. %d'%i 
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains)<max_strain_num:
            strain_mask = strain_index.get_indexer(strains)
            
            # Already normalized snps
            snps = data_g['norm_snps'][...]
            freqs = data_g['freqs'][...]
            mafs = sp.minimum(freqs,1-freqs)
            maf_mask = mafs>min_maf
            snps = snps[maf_mask]

            # Calculating the average genotype matrix (snps rows, individuals collumns)
            snps = average_genotype_matrix(snps)

            if len(snps)==0:
                continue
            K_snps_slice = K_snps[strain_mask]
            K_snps_slice[:,strain_mask] += sp.dot(snps.T,snps)
            K_snps[strain_mask] = K_snps_slice
            counts_mat_snps_slice = counts_mat_snps[strain_mask]
            counts_mat_snps_slice[:,strain_mask] += len(snps)
            counts_mat_snps[strain_mask] = counts_mat_snps_slice
              
    K_snps  = K_snps/counts_mat_snps  #element-wise division
    print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))
    correlation_plot(K_snps, wrt = False)
    #K_snps = pd.DataFrame(K_snps)
    #K_snps.to_csv('kinship_all_core.csv', header = ordered_strains) # The first collumn of the data contains the way strains are sorted
   
    tuple_index_kinship = (strain_index, K_snps)
    return(tuple_index_kinship)
kinship_all_genes()

def simple_mantel_nod_genes_nod_genes(max_strain_num=198,
                            maf=0.1,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Nod genes versus nod genes"""
    nod_genes = parse_nod()
    print nod_genes.keys()
    # Decoding the nod gene names
    nod_list = []
    for i in nod_genes.keys():
        nod_list.append(str(i).decode("utf-8"))

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()

    gene_grm_dict = {}
    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(nod_list), len(nod_list)))
    cor_matrix = pd.DataFrame(cor_matrix, index=nod_genes.values(), columns=nod_genes.values())

    for i in nod_list:
        try:
            data = h5f[i]
            total_snps_1 = data['snps'][...].T
            print 'The gene %s has %s snps' % (nod_genes[int(i)], total_snps_1.shape)
            filt = minor_allele_filter(total_snps_1, maf)
            print 'After filter MAF > %f has: %s' % (maf, filt.shape[1])
        except KeyError:
            'The nod gene %s is not in our subset of the data' % nod_genes[int(i)]

    for i, gg1 in enumerate(nod_list):
        try:
            strains_1 = h5f[gg1]['strains'][...]
        except KeyError:
            print 'The nod gene %s is not in our subset of the data' % nod_genes[int(gg1)]
            continue
        assert len(np.unique(strains_1)) == len(strains_1), 'There are multiple instances of the same strain in the data?'

        for j, gg2 in enumerate(nod_list):

            try:
                strains_2 = h5f[gg2]['strains'][...]
            except KeyError:
                # print 'The nod gene %s is not in our subset of the data' % nod_genes[int(gg2)]
                continue
            assert len(np.unique(strains_2)) == len(strains_2), 'There are multiple instances of the same strain in the data?'

            # This works only if we assume that strains_2 and strains_1 are ordered beforehand.  Are they? They are.
            strain_mask_1 = np.in1d(strains_1, strains_2, assume_unique=True)
            fitered_strains_1 = strains_1[strain_mask_1]
            strain_mask_2 = np.in1d(strains_2, fitered_strains_1, assume_unique=True)
            fitered_strains_2 = strains_2[strain_mask_2]

            # Finding the strains in common to the overall kinship (unpacking a tuple)
            (kinship_index, kinship) = kinship_all_genes()
            print len(strains_1)
            print len(kinship_index)
            strain_mask_kinship = np.in1d(kinship_index, strains_1, assume_unique=True)
            common_strains = strain_mask_kinship[strain_mask_kinship]
            kinship_subset = kinship[common_strains,:]
            kinship_subset = kinship_subset[:,common_strains]
            print kinship_subset.shape

            # 2. Calculate A, the cholesky decomp of the inverse of the GRM.
            #print("Finding inverse and sqrt of covariance matrix...")

            #inv_cov_sqrt = linalg.cholesky(linalg.inv(cov))

            #print inv_cov_sqrt

            # 3. Calculate the pseudo-SNPs (x*A)
            #print("Calculating pseudo SNPS...")
            #pseudo_snps = np.column_stack(np.dot(inv_cov_sqrt, col) for col in full_genotype_matrix.T)
            #del full_genotype_matrix


            # Only use the following code if you have all strains (or a fixed common set of strains).
#             if gg1 not in gene_grm_dict:
#                 data_g1 = h5f[gg1]
#                 total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns
#
#                 # Calculating GRM
#                 total_snps_1 = minor_allele_filter(total_snps_1, maf)
#                 grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
#                 gene_grm_dict[str(gg1)] = {'grm':grm_1}
#
#             if gg2 not in gene_grm_dict:
#                 data_g2 = h5f[gg2]
#                 total_snps_2 = data_g2['snps'][...].T  # strains in the rows, snps in the columns
#
#                 # Calculating GRM
#                 total_snps_2 = minor_allele_filter(total_snps_2, maf)
#                 grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])
#                 gene_grm_dict[str(gg2)] = {'grm':grm_2}

            data_g1 = h5f[gg1]
            total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns
            total_snps_1 = total_snps_1[strain_mask_1,:]  # Assuming it's number of SNPs x number of strains (Maria: I transposed so it is in the other way around)

            # Calculating GRM
            total_snps_1 = minor_allele_filter(total_snps_1, maf)
            grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
            gene_grm_dict[str(gg1)] = {'grm':grm_1}

            data_g2 = h5f[gg2]
            total_snps_2 = data_g2['snps'][...].T  # strains in the rows, snps in the columns
            total_snps_2 = total_snps_2[strain_mask_2,:]

            # Calculating GRM
            total_snps_2 = minor_allele_filter(total_snps_2, maf)
            grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])
            gene_grm_dict[str(gg2)] = {'grm':grm_2}

            # Calculating correlation and covariance based on the common subset of strains
            grm_1 = gene_grm_dict[str(gg1)]['grm']
            flat_grm_1 = grm_1.flatten()
            norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
            norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))

            grm_2 = gene_grm_dict[str(gg2)]['grm']
            flat_grm_2 = grm_2.flatten()
            norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
            norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))

            # Built in function, it returns correlation coefficient and the p-value for testing non-correlation
            r = pearsonr(norm_flat_grm1, norm_flat_grm2)
            cor_matrix[nod_genes[int(gg1)]][nod_genes[int(gg2)]] = r[0]

            # Calculating against the overall kinship
            print gg1
            print pearsonr(norm_flat_grm1, kinship_subset.flatten())

    correlation_plot(cor_matrix)
    cor_matrix.to_csv('Mantel_test_nod_all_maf_1.csv', header=True)
    #correlation_plot(cor_matrix)

#simple_mantel_nod_genes_nod_genes()
