'''Implement a simpler version of the Mantel test, ignoring population structure.'''

# In this script we calculate the kinship matrix for each core gene and do a simple mantel test (correlation of matrix) between pairs of genes

# MantelTest v1.2.10
# http://jwcarr.github.io/MantelTest/
#
# Copyright (c) 2014-2016 Jon W. Carr
# Licensed under the terms of the MIT License

import numpy as np
from itertools import permutations
from scipy import spatial, stats
import scipy as sp
import h5py
import pandas as pd
import pylab as pl
from scipy.stats.stats import pearsonr 


def minor_allele_filter(gene_matrix, maf):
    '''Filtering for minor allele frequency, it assumes that the matrix is binary and that is n x m, where columns are markers (m).'''

    freqs = sp.mean(gene_matrix, 0) 
    mafs = sp.minimum(freqs, 1 - freqs)

    maf_filter = mafs > maf
    matrix_mafs = gene_matrix[:, maf_filter]

    # Normalizing the columns:
    norm_matrix = (matrix_mafs - np.mean(matrix_mafs, axis=0)) / np.std(matrix_mafs, axis=0)
    return(norm_matrix)

def corr2_coeff(A,B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);

    # Finally get corr coeff
    return np.dot(A_mA,B_mB.T)/np.sqrt(np.dot(ssA[:,None],ssB[None]))

def simple_intergenic_ld_core(max_strain_num=198,
                            maf=0.2,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    
    core_genes = []
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(set(strains)) == max_strain_num:
            core_genes.append(gg)
    
    r_scores = []
    gene_names = []
    gene_grm_dict = {}

    for i, gg1 in enumerate(core_genes):
        data_g1 = h5f[gg1]
        total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns 
        #print total_snps_1.shape
        total_snps_1 = minor_allele_filter(total_snps_1, 0.1)

        ''' 3. Calculate the Kinship matrix for each gene '''
        grm = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])

        flat_grm = grm.flatten()
        norm_flat_grm1 = flat_grm - flat_grm.mean() / sp.sqrt(sp.dot(flat_grm, flat_grm))
        #print np.var(norm_flat_grm1)
        #print np.mean(norm_flat_grm1)
        gene_grm_dict[gg1] = {'grm':grm , 'norm_flat_grm':norm_flat_grm1}
        
        for j, gg2 in enumerate(core_genes):
            if i > j:
                            
                norm_flat_grm2 = gene_grm_dict[gg2]['norm_flat_grm']
                covariance = sp.dot(norm_flat_grm1, norm_flat_grm2) 
                var1 = np.sum(abs(norm_flat_grm1 - norm_flat_grm1.mean())**2)
                var2 = np.sum(abs(norm_flat_grm2 - norm_flat_grm2.mean())**2)
                r = covariance/sp.sqrt(var1 * var2)
                
                # Checking the values with a scipy built function:
                #r_bel = pearsonr(norm_flat_grm1, norm_flat_grm2)
                #print round(r, 5) == round(r_bel[0], 5)

                r_scores.append(r)
                gene_names.append(gg1 + gg2)

    LD_stats = pd.DataFrame(
    {'r_scores': r_scores,
    'gene_names': gene_names})

    LD_stats.to_csv('test.csv', header=True)

    return LD_stats

simple_intergenic_ld_core()


def simple_intergenic_ld_nod_genes(max_strain_num=198,
                            maf=0.2,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    nod_genes = {'nodA': 8048, 'nodB': 9911, 'nodC': 10421, 'nodD':7218, 'nodE': 4140, 'nodF': 4139, 'nifA': 4128, 'nifB': 4129, 'nifH': 4123, 'nifK': 4121, 'fixA': 4124, 'fixB': 4125, 'fixC': 4126, 'fixN': 4130}
    
    nod_list = []
    for i in nod_genes.values():
      nod_list.append(str(i).decode("utf-8"))

    print nod_genes.values()
    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    
    core_genes = []
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(set(strains)) == max_strain_num:
            core_genes.append(gg)
    
    r_scores = []
    gene_names = []
    gene_grm_dict = {}

    for i, gg1 in enumerate(nod_list):
        data_g1 = h5f[gg1]
        total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns 
        total_snps_1 = minor_allele_filter(total_snps_1, 0.1)

        #Calculate the Kinship matrix for each gene
        grm = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])

        flat_grm = grm.flatten()
        norm_flat_grm1 = flat_grm - flat_grm.mean() / sp.sqrt(sp.dot(flat_grm, flat_grm))
        gene_grm_dict[gg1] = {'grm':grm , 'norm_flat_grm':norm_flat_grm1}
        
        for j, gg2 in enumerate(core_genes):
            if i > j:
                            
                norm_flat_grm2 = gene_grm_dict[gg2]['norm_flat_grm']
                covariance = sp.dot(norm_flat_grm1, norm_flat_grm2) 
                var1 = np.sum(abs(norm_flat_grm1 - norm_flat_grm1.mean())**2)
                var2 = np.sum(abs(norm_flat_grm2 - norm_flat_grm2.mean())**2)
                r = covariance/sp.sqrt(var1 * var2)
                
                # Checking the values with a scipy built in function:
                #r_bel = pearsonr(norm_flat_grm1, norm_flat_grm2)
                #print round(r, 5) == round(r_bel[0], 5)

                r_scores.append(r)
                gene_names.append(gg1 + gg2)

    LD_stats = pd.DataFrame(
    {'r_scores': r_scores,
    'gene_names': gene_names})

    LD_stats.to_csv('test.csv', header=True)

    return LD_stats

#simple_intergenic_ld_nod_genes()


def simple_intergenic_ld_accessory(specific_genes,
                            max_strain_num=198,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    # """Gives a specific list of genes and calculate LD of these genes with all"""


    nod_names = open(specific_genes)

    specific_genes = [line.rstrip('\n') for line in nod_names]
    print specific_genes

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    all_strains = set()
    
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(set(strains)) == max_strain_num:
            all_strains = set(strains).union(all_strains)
    num_strains = len(all_strains)
    print 'Found %d "distinct" strains' % num_strains
    
    ordered_strains = sorted(list(all_strains))
    strain_index = pd.Index(ordered_strains)

    count = 0
    r_scores = []
    p_values = []
    z_scores = []
    gene_names = []

    ag = h5f['alignments']
    print ag['8048']

    for i, gg1 in enumerate(gene_groups):
        for j, gg1 in enumerate(gene_groups):
            if i > j:
            
                data_g1 = h5f[gg1]
                data_g2 = h5f[gg2]  # tuple

                strains_1 = data_g1['strains'][...]
                strains_2 = data_g2['strains'][...]
                            
                # Indexes of the strains in the big GRM matrix
                strain_mask_1 = strain_index.get_indexer(strains_1)
                strain_mask_2 = strain_index.get_indexer(strains_2)
                
                set_1, set_2 = set(strain_mask_1), set(strain_mask_2)
                intersec = list(set_1 & set_2)

                # Indexes of the strains in their own matrix
                    
                # list1 determines the ordering
                olist1 = [i for i, item in enumerate(strains_1) if item in set(strains_2)]

                # list2 determines the ordering
                olist2 = [i for i, item in enumerate(strains_2) if item in set(strains_1)]

                # Take the subset of shared snps of each gene
                total_snps_1 = data_g1['norm_snps'][...].T  # strains in the rows, snps in the collumns 
                common_snps_1 = total_snps_1[olist1, :]

                total_snps_2 = data_g2['norm_snps'][...].T
                common_snps_2 = total_snps_2[olist2, :]

                ''' 3. Calculate the Kinship matrix for each gene '''

                pseudo_snps_1 = np.dot(common_snps_1, common_snps_1) / common_snps_1[1]
                pseudo_snps_2 = np.dot(common_snps_2, common_snps_2) / common_snps_2[1]

                r, p, z = Mantel.mantel_test(grm_1, grm_2, perms=1, method='spearman', tail='two-tail')

                # print(r,p,z)
                r_scores.append(r)
                p_values.append(p)
                z_scores.append(z)
                gene_names.append(gg1[0][0:5] + gg2)

        LD_stats = pd.DataFrame(
        {'r_scores': r_scores,
        'p_values': p_values,
        'z_scores': z_scores,
        'gene_names': gene_names})

        LD_stats.to_csv('test.csv', header=True)


    return LD_stats

