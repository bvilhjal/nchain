'''Implement a simpler version of the Mantel test, ignoring population structure.'''

# In this script we calculate the kinship matrix for each core gene and do a simple mantel test (correlation of matrix) between pairs of genes

# MantelTest v1.2.10
# http://jwcarr.github.io/MantelTest/
#
# Copyright (c) 2014-2016 Jon W. Carr
# Licensed under the terms of the MIT License

import numpy as np
from scipy import spatial, stats
import scipy as sp
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr 
import seaborn as sns
import time

def minor_allele_filter(gene_matrix, maf):
    '''Filtering for minor allele frequency, it assumes that the matrix is binary and that is n x m, where columns are markers (m).'''

    freqs = sp.mean(gene_matrix, 0) 
    mafs = sp.minimum(freqs, 1 - freqs)

    maf_filter = mafs > maf
    matrix_mafs = gene_matrix[:, maf_filter]

    # Normalizing the columns:
    norm_matrix = (matrix_mafs - np.mean(matrix_mafs, axis=0)) / np.std(matrix_mafs, axis=0)
    return(norm_matrix)
   
def correlation_plot(df):
  corr = df.corr()
  mask = np.zeros_like(corr)
  #mask[np.triu_indices_from(mask)] = True

  # Set up the matplotlib figure
  f, ax = plt.subplots(figsize=(12, 9))

  # Clustering with seaborn
  with sns.axes_style("white"):
    ax = sns.heatmap(corr, mask=mask, vmax=.3, square=True, annot=True, annot_kws={"size": 9})
  f.tight_layout()
  plt.show()

def simple_intergenic_ld_core(max_strain_num=198,
                            maf=0.2,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    t0 = time.time()
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

    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(core_genes[0:100]), len(core_genes[0:100])))
    cor_matrix = pd.DataFrame(cor_matrix, index = core_genes[0:100], columns = core_genes[0:100])

    for i, gg1 in enumerate(core_genes[0:100]):
        data_g1 = h5f[gg1]
        total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns 
        total_snps_1 = minor_allele_filter(total_snps_1, 0.1)

        ''' 3. Calculate the Kinship matrix for each gene '''
        grm = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])

        flat_grm = grm.flatten()
        norm_flat_grm1 = flat_grm - flat_grm.mean() / sp.sqrt(sp.dot(flat_grm, flat_grm))
        gene_grm_dict[gg1] = {'grm':grm , 'norm_flat_grm':norm_flat_grm1}
        
        for j, gg2 in enumerate(core_genes[0:1000]):
            if i > j:
                            
                norm_flat_grm2 = gene_grm_dict[gg2]['norm_flat_grm']
                covariance = sp.dot(norm_flat_grm1, norm_flat_grm2) 
                var1 = np.sum(abs(norm_flat_grm1 - norm_flat_grm1.mean())**2)
                var2 = np.sum(abs(norm_flat_grm2 - norm_flat_grm2.mean())**2)
                r = covariance/sp.sqrt(var1 * var2)

                cor_matrix[gg2][gg1] = r
                
                # Checking the values with a scipy built function:
                #r_bel = pearsonr(norm_flat_grm1, norm_flat_grm2)
                #print round(r, 5) == round(r_bel[0], 5)

    print cor_matrix
    LD_stats.to_csv('test.csv', header=True)
    t1 = time.time()

    total = t1-t0
    print 'total amount of time consumed is %f' % total
    return LD_stats

simple_intergenic_ld_core()


def simple_intergenic_ld_nod_genes(max_strain_num=198,
                            maf=0.2,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    nod_genes = {8048:'nodA', 9911: 'nodB', 10421: 'nodC', 7218: 'nodD', 4140: 'nodE', 4139: 'nodF', 10588 : 'nodI', 4134: 'nodJ', 4141: 'nodL', 4142: 'nodM', 4144: 'nodX', 4143: 'nodN', 4128: 'nifA', 4129: 'nifB', 4122: 'nifD', 4123: 'nifH', 4121: 'nifK', 4119: 'nifN', 4124: 'fixA', 4125: 'fixB', 4126:'fixC', 4130:'fixN', 4127: 'fixX'}
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
    
    r_scores = []
    gene_names = []
    gene_grm_dict = {}

    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(nod_list), len(core_genes[0:100])))
    cor_matrix = pd.DataFrame(cor_matrix, index = nod_genes.values(), columns = core_genes[0:100])

    for i, gg1 in enumerate(nod_list):
        try:
          strains_1 = h5f[gg1]['strains'][...]
        except KeyError:
          print 'The nod gene %s is not in our subset of the data' % nod_genes[int(gg1)]
          continue
        
        for j, gg2 in enumerate(core_genes[0:100]):

          strains_2 = h5f[gg2]['strains'][...]
        
          set_1, set_2 = set(strains_1), set(strains_2)
          intersec = list(set_1 & set_2)

          strain_mask_2 = []
          strain_mask_1 = []

          for i in intersec:
              strain_mask_2.append(np.unique(strains_2).tolist().index(i))
              strain_mask_1.append(np.unique(strains_1).tolist().index(i))

          strain_mask_2 = sorted(strain_mask_2)
          strain_mask_1 = sorted(strain_mask_1)

          if gg1 not in gene_grm_dict:
            data_g1 = h5f[gg1]
            total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns 

            # Calculating GRM 
            total_snps_1 = minor_allele_filter(total_snps_1, 0.1)
            grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
            gene_grm_dict[str(gg1)] = {'grm':grm_1}

          if gg2 not in gene_grm_dict:
            data_g2 = h5f[gg2]
            total_snps_2 = data_g2['snps'][...].T  # strains in the rows, snps in the columns 
          
            # Calculating GRM 
            total_snps_2 = minor_allele_filter(total_snps_2, 0.1)
            grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])
            gene_grm_dict[str(gg2)] = {'grm':grm_2}

          # Calculating correlation and covariance based on the common subset of strains
          grm_1 = gene_grm_dict[str(gg1)]['grm']
          sub_grm_1 = grm_1[strain_mask_1,strain_mask_1]
          flat_grm_1 = sub_grm_1.flatten()
          norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean() / sp.sqrt(sp.dot(flat_grm_1, flat_grm_1))
          
          grm_2 = gene_grm_dict[str(gg2)]['grm']
          sub_grm_2 = grm_2[strain_mask_2,strain_mask_2]
          flat_grm_2 = sub_grm_2.flatten()
          norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean() / sp.sqrt(sp.dot(flat_grm_2, flat_grm_2))

          # Built in function, it returns correlation coefficient and the p-value for testing non-correlation
          r_bel = pearsonr(norm_flat_grm1, norm_flat_grm2)

          # Indexing column (gg2) and row (gg1)
          cor_matrix[gg2][nod_genes[int(gg1)]] += r_bel[0]
          #covariance = sp.dot(norm_flat_grm1, norm_flat_grm2) 
          #var1 = np.sum(abs(norm_flat_grm1 - norm_flat_grm1.mean())**2)
          #var2 = np.sum(abs(norm_flat_grm2 - norm_flat_grm2.mean())**2)
          #r = covariance/sp.sqrt(var1 * var2)
 
          r_scores.append(r_bel[0])
          gene_names.append(nod_genes[int(gg1)] +'_'+ gg2)          

    cor_matrix.to_csv('Mantel_test_nod_all.csv', header = True)     
    return cor_matrix

#simple_intergenic_ld_nod_genes()

def simple_mantel_nod_genes_nod_genes(max_strain_num=198,
                            maf=0.2,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    #nod_genes = {8048:'nodA', 9911: 'nodB', 10421: 'nodC', 7218: 'nodD', 4140: 'nodE', 4139: 'nodF', 10588 : 'nodI', 4134: 'nodJ', 4141: 'nodL', 4142: 'nodM', 4144: 'nodX', 4143: 'nodN', 4128: 'nifA', 4129: 'nifB', 4122: 'nifD', 4123: 'nifH', 4121: 'nifK', 4119: 'nifN', 4124: 'fixA', 4125: 'fixB', 4126:'fixC', 4130:'fixN', 4127: 'fixX'}
    nod_genes = {4140: 'nodE', 4139: 'nodF', 4134: 'nodJ', 4141: 'nodL', 4142: 'nodM', 4144: 'nodX', 4143: 'nodN', 4128: 'nifA', 4129: 'nifB', 4122: 'nifD', 4123: 'nifH', 4121: 'nifK', 4119: 'nifN', 4124: 'fixA', 4125: 'fixB', 4126:'fixC', 4127: 'fixX'}
    # Decoding the nod gene names
    nod_list = []
    for i in nod_genes.keys():
      nod_list.append(str(i).decode("utf-8"))
    #print nod_genes.values()

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    
    r_scores = []
    gene_names = []
    gene_grm_dict = {}

    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(nod_list), len(nod_list)))
    cor_matrix = pd.DataFrame(cor_matrix, index = nod_genes.values(), columns = nod_genes.values())

    for i in nod_list:
      try:
          data = h5f[i]
          total_snps_1 = data['snps'][...].T 
          print 'The gene %s has %s snps' % (nod_genes[int(i)], total_snps_1.shape)
          filt = minor_allele_filter(total_snps_1, 0.1)
          print 'After filer MAF > 0.1 has: %s' %  (filt.shape[1])
      except KeyError:
          'The nod gene %s is not in our subset of the data' % nod_genes[int(i)]


    for i, gg1 in enumerate(nod_list):
        try:
          strains_1 = h5f[gg1]['strains'][...]
        except KeyError:
          #print 'The nod gene %s is not in our subset of the data' % nod_genes[int(gg1)]
          continue
        
        for j, gg2 in enumerate(nod_list):

          try:
            strains_2 = h5f[gg2]['strains'][...]
          except KeyError:
           #print 'The nod gene %s is not in our subset of the data' % nod_genes[int(gg2)]
            continue
        
          set_1, set_2 = set(strains_1), set(strains_2)
          intersec = list(set_1 & set_2)

          strain_mask_2 = []
          strain_mask_1 = []

          for i in intersec:
              strain_mask_2.append(np.unique(strains_2).tolist().index(i))
              strain_mask_1.append(np.unique(strains_1).tolist().index(i))

          strain_mask_2 = sorted(strain_mask_2)
          strain_mask_1 = sorted(strain_mask_1)

          if gg1 not in gene_grm_dict:
            data_g1 = h5f[gg1]
            total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns 

            # Calculating GRM 
            total_snps_1 = minor_allele_filter(total_snps_1, 0.2)
            grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
            gene_grm_dict[str(gg1)] = {'grm':grm_1}

          if gg2 not in gene_grm_dict:
            data_g2 = h5f[gg2]
            total_snps_2 = data_g2['snps'][...].T  # strains in the rows, snps in the columns 
          
            # Calculating GRM 
            total_snps_2 = minor_allele_filter(total_snps_2, 0.2)
            grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])
            gene_grm_dict[str(gg2)] = {'grm':grm_2}

          # Calculating correlation and covariance based on the common subset of strains
          grm_1 = gene_grm_dict[str(gg1)]['grm']
          sub_grm_1 = grm_1[strain_mask_1,strain_mask_1]
          flat_grm_1 = sub_grm_1.flatten()
          norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean() / sp.sqrt(sp.dot(flat_grm_1, flat_grm_1))
          
          grm_2 = gene_grm_dict[str(gg2)]['grm']
          sub_grm_2 = grm_2[strain_mask_2,strain_mask_2]
          flat_grm_2 = sub_grm_2.flatten()
          norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean() / sp.sqrt(sp.dot(flat_grm_2, flat_grm_2))

          # Built in function, it returns correlation coefficient and the p-value for testing non-correlation
          r_bel = pearsonr(norm_flat_grm1, norm_flat_grm2)

          cor_matrix[nod_genes[int(gg1)]][nod_genes[int(gg2)]] += r_bel[0]

          #covariance = sp.dot(norm_flat_grm1, norm_flat_grm2) 
          #var1 = np.sum(abs(norm_flat_grm1 - norm_flat_grm1.mean())**2)
          #var2 = np.sum(abs(norm_flat_grm2 - norm_flat_grm2.mean())**2)
          #r = covariance/sp.sqrt(var1 * var2)
          
          #print r    
          r_scores.append(r_bel[0])
          #gene_names.append(nod_genes[int(gg1)] +'_'+ nod_genes[int(gg2)])   


    cor_matrix.to_csv('Mantel_test_nod_all_maf_2.csv', header = True)
    correlation_plot(cor_matrix)
    #LD_stats = pd.DataFrame(
    #{'r_scores': r_scores,
    #'gene_names': gene_names})

    #LD_stats.to_csv('test_nod_genes.csv', header=True)

    #return LD_stats

#simple_mantel_nod_genes_nod_genes()