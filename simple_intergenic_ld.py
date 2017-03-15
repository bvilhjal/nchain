'''Implement a simpler version of the Mantel test, ignoring population structure.'''

# In this script we calculate the kinship matrix for each core gene and do a simple mantel test (correlation of matrix) between pairs of genes

import numpy as np
import scipy as sp
import h5py
import bottleneck as bn
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import seaborn as sns
import time
from collections import OrderedDict
from numpy import linalg
import pylab as pl
import glob 

def parse_nod():
    nod_genes = OrderedDict()
    # Full matrix:
    # nod_genes = {4144: 'nodX', 4143: 'nodN', 4142: 'nodM', 4141: 'nodL',10151: 'nodR', 4140: 'nodE', 4139: 'nodF', 4138: 'nodD', 4137: 'nodA', 8237: 'nodB', 4136: 'nodC', 4135: 'nodI', 4134: 'nodJ', 4131: 'fixT', 4130:'fixN', 4129: 'nifB', 4128: 'nifA', 4127: 'fixX', 4126:'fixC', 4125: 'fixB', 4124: 'fixA', 4123: 'nifH', 4122: 'nifD', 4121: 'nifK', 4120: 'nifE', 4119: 'nifN', 2448: 'rpoB', 2140: 'recA'}
    nod_genes = OrderedDict([(4144, 'nodX'), (4143, 'nodN'), (4142, 'nodM'), (4141, 'nodL'), (4140, 'nodE'), (4139, 'nodF'), (4138, 'nodD'), (4137, 'nodA'), (4136, 'nodC'), (4135, 'nodI'), (4134, 'nodJ'), (4129, 'nifB'), (4128, 'nifA'), (4127, 'fixX'), (4126, 'fixC'), (4125, 'fixB'), (4124, 'fixA'), (4123, 'nifH'), (4122, 'nifD'), (4121, 'nifK'), (4120, 'nifE'), (2448, 'rpoB'), (2140, 'recA')])
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
    '''Computes the average genotype matrix, it assumes that the input matrix (Markers/M in cols and individuals/N in rows)'''

    # Normalizing the rows:
    average_X = (X - np.mean(X, axis=1))
    return(average_X)

def correlation_plot(df, wrt=True):
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(12, 9))

    # Clustering with seaborn
    with sns.axes_style("white"):
        ax = sns.heatmap(df, square=True, annot=wrt, annot_kws={"size": 9}, cmap="RdYlGn", vmin=0, vmax=1)
    f.tight_layout()
    plt.show()

def kinship_all_genes(snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5',
                 min_maf=0.05,
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
        if len(strains) < max_strain_num:
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
    correlation_plot(K_snps, wrt=False)

    tuple_index_kinship = (strain_index, K_snps)
    return(tuple_index_kinship)
# kinship_all_genes()

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

    # Opening the kinship matrix
    # (kinship_index, kinship) = kinship_all_genes()

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
            total_snps_1 = total_snps_1[strain_mask_1, :]

            total_snps_1 = minor_allele_filter(total_snps_1, maf)
            grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
            gene_grm_dict[str(gg1)] = {'grm':grm_1}

            data_g2 = h5f[gg2]
            total_snps_2 = data_g2['snps'][...].T  # strains in the rows, snps in the columns
            total_snps_2 = total_snps_2[strain_mask_2, :]

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
            # print pearsonr(norm_flat_grm1, kinship_subset.flatten())

    correlation_plot(cor_matrix)
    cor_matrix.to_csv('Mantel_test_nod_all_maf_1.csv', header=True)
    # correlation_plot(cor_matrix)

#simple_mantel_nod_genes_nod_genes()

def mantel_corrected_nod_genes(in_glob = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/'):
    """Take the structured corrected files and calculate mantel test for the nod genes"""
    parse_nod_genes = parse_nod() 


    #'group1000.npz', 'group2000.npz', 'group10.npz'
    nod_genes = ['group4144.npz', 'group4143.npz', 'group4142.npz', 'group4141.npz', 'group4140.npz', 'group4139.npz', 'group4138.npz', 'group4137.npz', 'group4136.npz', 'group4135.npz', 'group4134.npz', 'group4129.npz', 'group4128.npz', 'group4127.npz', 'group4126.npz', 'group4125.npz', 'group4124.npz', 'group4123.npz', 'group4122.npz', 'group4121.npz', 'group4120.npz', 'group2448.npz', 'group2140.npz']
    
    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(nod_genes), len(nod_genes)))
    cor_matrix = pd.DataFrame(cor_matrix, index=parse_nod_genes.values(), columns=parse_nod_genes.values())

    #This is the gene SNPs matrix
    genes = []
    for f in nod_genes:
        with np.load(in_glob + f) as data:
            # Creating a tuple
            genes.append((f, data["matrix"], data["strains"])) 

    for gg1 in genes:
        for gg2 in genes:

            gene1, total_snps_1, strains_1 = (gg1)
            gene2, total_snps_2, strains_2 = (gg2)

            # This works only if we assume that strains_2 and strains_1 are ordered beforehand.  Are they? They are.
            strain_mask_1 = np.in1d(strains_1, strains_2, assume_unique=True)
            fitered_strains_1 = strains_1[strain_mask_1]
            strain_mask_2 = np.in1d(strains_2, fitered_strains_1, assume_unique=True)
            fitered_strains_2 = strains_2[strain_mask_2]

            total_snps_1 = total_snps_1[strain_mask_1, :]
            grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])

            total_snps_2 = total_snps_2[strain_mask_2, :]
            grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])

            # Calculating correlation and covariance based on the common subset of strains
            flat_grm_1 = grm_1.flatten()
            norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
            norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))

            flat_grm_2 = grm_2.flatten()
            norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
            norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))

            # Built in function, it returns correlation coefficient and the p-value for testing non-correlation
            print parse_nod_genes[int(gene1[5:9])] +'_'+ parse_nod_genes[int(gene2[5:9])]
            r = pearsonr(norm_flat_grm1, norm_flat_grm2)
            r2 = pearsonr(flat_grm_1, flat_grm_2)
            print r[0]
            print r2[0]
            cor_matrix[parse_nod_genes[int(gene1[5:9])]][parse_nod_genes[int(gene2[5:9])]] = r[0]
    correlation_plot(cor_matrix)

mantel_corrected_nod_genes()

def simple_intergenic_ld_nod_genes(max_strain_num=198,
                            maf=0.1,
                            amount = 10,
                            snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    nod_genes = parse_nod()
    # Decoding the nod gene names
    nod_list = []
    for i in nod_genes.keys():
        nod_list.append(str(i).decode("utf-8"))
    print nod_genes.values()

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()

    print h5f[nod_list[0]]
    core_genes = []
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(set(strains)) == max_strain_num:
            core_genes.append(gg)
    core_genes = sorted(core_genes)
    gene_grm_dict = {}

    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(nod_list), len(core_genes[0::])))
    cor_matrix = pd.DataFrame(cor_matrix, index=nod_genes.values(), columns=core_genes[0::])

    for gg1 in nod_list:
        try:
            strains_1 = h5f[gg1]['strains'][...]
        except KeyError:
            print 'The nod gene %s is not in our subset of the data' % nod_genes[int(gg1)]
            continue

        for gg2 in core_genes[0::]:

            strains_2 = h5f[gg2]['strains'][...]

            # This works only if we assume that strains_2 and strains_1 are ordered beforehand.  Are they? They are.
            strain_mask_1 = np.in1d(strains_1, strains_2, assume_unique=True)
            fitered_strains_1 = strains_1[strain_mask_1]
            strain_mask_2 = np.in1d(strains_2, fitered_strains_1, assume_unique=True)
            fitered_strains_2 = strains_2[strain_mask_2]

            data_g1 = h5f[gg1]
            total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns
            total_snps_1 = total_snps_1[strain_mask_1, :]

            total_snps_1 = minor_allele_filter(total_snps_1, maf)
            grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
            gene_grm_dict[str(gg1)] = {'grm':grm_1}

            data_g2 = h5f[gg2]
            total_snps_2 = data_g2['snps'][...].T  # strains in the rows, snps in the columns
            total_snps_2 = total_snps_2[strain_mask_2, :]

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
            cor_matrix[gg2][nod_genes[int(gg1)]] += r[0]
            #print r
    
    #correlation_plot(cor_matrix, wrt = false)
    cor_matrix.to_csv('Mantel_test_nod_all_core.csv', header=True)
    return cor_matrix

#simple_intergenic_ld_nod_genes()