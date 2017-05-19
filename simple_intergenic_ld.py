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
from sys import argv

# Setting the directory of the data depending on what computer I am  

if argv[1] != 'mac':
    snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'
    out_dir='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps_test/'
    figure_dir='C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Figures/'
    in_glob = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps_test/'

else:
    snps_file='/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5'
    out_dir='/Users/PM/Desktop/PHD_incomplete/Methods/Intergenic_LD/corrected_snps_test/'
    figure_dir='/Users/PM/Desktop/PHD_incomplete/nchain/Figures/'  
    in_glob = '/Users/PM/Desktop/PHD_incomplete/Methods/Intergenic_LD/corrected_snps_test/'

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

def correlation_plot(df, wrt=True, fig_name = 'mantel_test.png', show = False,
                    figure_dir= figure_dir):
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(12, 9))

    # Clustering with seaborn
    with sns.axes_style("white"):
        ax = sns.heatmap(df, square=True, annot=wrt, annot_kws={"size": 9}, cmap="RdYlGn", vmin=0, vmax=1)
    f.tight_layout()
    plt.savefig(figure_dir + fig_name)
    if show == True:
        plt.show()

def simple_mantel_nod_genes_nod_genes(max_strain_num=198,
                            maf=0.1,
                           snps_file= snps_file):
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

    #correlation_plot(cor_matrix)
    #cor_matrix.to_csv('Mantel_test_nod_all_maf_1.csv', header=True)
    return(cor_matrix)
    # correlation_plot(cor_matrix)

#simple_mantel_nod_genes_nod_genes()

def mantel_corrected_nod_genes(in_glob = in_glob,
                                min_maf = 0.1,
                                snps_file = snps_file, 
                                fig_name = 'default.pdf',
                                figure_dir= figure_dir,
                                slicing = True):
    """Take the structured corrected files and calculate mantel test for the nod genes"""
    parse_nod_genes = parse_nod() 

    nod_genes = ['group4144.npz', 'group4143.npz', 'group4142.npz', 'group4141.npz', 'group4140.npz', 'group4139.npz', 'group4138.npz', 'group4137.npz', 'group4136.npz', 'group4135.npz', 'group4134.npz', 'group4129.npz', 'group4128.npz', 'group4127.npz', 'group4126.npz', 'group4125.npz', 'group4124.npz', 'group4123.npz', 'group4122.npz', 'group4121.npz', 'group4120.npz', 'group2448.npz', 'group2140.npz']
    
    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(nod_genes), len(nod_genes)))
    cor_matrix = pd.DataFrame(cor_matrix, index=parse_nod_genes.values(), columns=parse_nod_genes.values())

    #This is the gene SNPs matrix
    genes = []
    for f in nod_genes:
        with np.load(in_glob + f) as data:
            # Creating a tuple of gene name, followed by snp matrix, strains mask, maf mask
            genes.append((f, data["matrix"], data["strains"])) 

    for gg1 in genes:
        for gg2 in genes:

            gene1, total_snps_1, strains_1 = (gg1)
            gene2, total_snps_2, strains_2 = (gg2)

            if slicing:
            
                # This works only if we assume that strains_2 and strains_1 are ordered beforehand.  Are they? They are.
                strain_mask_1 = np.in1d(strains_1, strains_2, assume_unique=True)
                fitered_strains_1 = strains_1[strain_mask_1]
                strain_mask_2 = np.in1d(strains_2, fitered_strains_1, assume_unique=True)

                # Shuffling randomly
                #strain_mask_1 = random.shuffle(strain_mask_1)
                #strain_mask_2 = random.shuffle(strain_mask_2)

                total_snps_1 = total_snps_1[strain_mask_1, :]

                grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])

                total_snps_2 = total_snps_2[strain_mask_2, :]
                grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])

            else:
                grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
                grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])

            # Calculating correlation and covariance based on the common subset of strains
            flat_grm_1 = grm_1.flatten()
            flat_grm_2 = grm_2.flatten()


            # Built in function, it returns correlation coefficient and the p-value for testing non-correlation
            r = pearsonr(flat_grm_1, flat_grm_2)

            cor_matrix[parse_nod_genes[int(gene1[5:9])]][parse_nod_genes[int(gene2[5:9])]] = r[0]
    correlation_plot(cor_matrix, show = False, fig_name = fig_name)
    return(cor_matrix)

#mantel_corrected_nod_genes()

def figure_comparison(corrected = 0, incorrect = 0):
    incorrected = simple_mantel_nod_genes_nod_genes()

    corrected = mantel_corrected_nod_genes()
    u = np.triu(corrected)
    l = np.tril(incorrected,k = -1)

    mix = u + l
    mix = pd.DataFrame(data = mix, columns = list(corrected), index = list(corrected))

    correlation_plot(mix)
#figure_comparison()

def robusteness_maf_simple():
    mafs = np.arange(0,0.5,0.1)
    for maf in mafs:

        cor_matrix = simple_mantel_nod_genes_nod_genes(maf = maf)
        name = 'mantel_test_wo_correction' +  str(maf) + '.pdf'
        correlation_plot(cor_matrix, fig_name = name)

#robusteness_maf_simple()

def simple_intergenic_ld_nod_genes(max_strain_num=100,
                            maf=0.1,
                            amount = 10,
                            snps_file= snps_file,
                            figure_dir= figure_dir):

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
        if len(set(strains)) > max_strain_num:
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
    
    #correlation_plot(cor_matrix, wrt = false)
    cor_matrix.to_csv('Mantel_test_nod_all_core.csv', header=True)
    return(cor_matrix)

#simple_intergenic_ld_nod_genes()