# 1 Chose Pacbio Data
# In the folder polished assemblies we have the needed files of each PacBio data to make the blast against the 5000 genes
from sys import argv 
from Bio.Blast.Applications import NcbiblastnCommandline
import h5py
import glob
import pandas as pd
import numpy as np
import scipy as sp
import re
from simple_intergenic_ld import correlation_plot
from sys import argv 
from itertools import izip
from collections import OrderedDict
from scipy.stats.stats import pearsonr
import os


if argv[1] == 'mac':
    snps_file='/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5'
    out_dir='/Users/PM/Desktop/PHD_incomplete/Methods/Intergenic_LD/corrected_snps_test/'
    meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'


def kinship_all_genes(snps_file= snps_file,
                 min_maf=0.10,
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

    tuple_index_kinship = (K_snps, ordered_strains)

    headers = list()
    maps = parse_pop_map()
    for i in ordered_strains:
        headers.append( 'SM' +maps[i]['sara_id'])

    np.savetxt('strains_order.csv', headers,  delimiter=",", fmt='%s')
    np.savetxt('kinship_maf_01.csv', K_snps, fmt='%.18e', delimiter=',', header = ','.join(headers), comments="")
    # Saving as a compressed numpy file:
    #np.savez_compressed("{}/{}".format(kinship_matrix), matrix=snps, strains=strains, maf = maf) 

    return(tuple_index_kinship)

def plot_dirty_PCA( figure_fn = 'pca.png', k_figure_fn = 'kinship_heatmap.png', title=None,
                   figure_dir = '/project/NChain/faststorage/rhizobium/ld/figures',strains=None):
    from scipy import linalg
    
    kinship_mat, strains = (kinship_all_genes())

    evals, evecs = linalg.eig(kinship_mat)  #PCA via eigen decomp
    evals[evals<0]=0
    sort_indices = sp.argsort(evals,)
    ordered_evals = evals[sort_indices]
    print ordered_evals[-10:]/sp.sum(ordered_evals)
    pc1,pc2 = evecs[:,sort_indices[-1]],evecs[:,sort_indices[-2]]
    pl.clf()
    
    if strains is not None:    
        ct_marker_map = {'DK':'*','UK':'^', 'F':'o', 'NA': 's'}
        gs_color_map = {'gsA':'#386CB0','gsB':'#FB8072', 'gsC':'#1B9E77', 'gsE': '#F0027F', 'gsD': '#984EA3'}
        pop_map = parse_pop_map()

        print pop_map
        for i, strain in enumerate(strains):
            print i
            d = pop_map.get(strain,'NA')
            if d=='NA':
                gs = 'NA'
                country = 'NA'
            else:
                gs = d['genospecies']
                country = d['country']
            pl.scatter(pc1[i],pc2[i], marker=ct_marker_map[country], c=gs_color_map[gs], alpha=0.3, s=100, edgecolor='none')
        for gs in gs_color_map:
            pl.scatter([], [], color=gs_color_map[gs], marker = 's', label=gs, s=100, edgecolor='none')
        for country in ct_marker_map:
            if country !='NA':
                pl.scatter([], [], color='k', marker = ct_marker_map[country], label=country, s=100, facecolors='none')

        
        pl.legend(scatterpoints=1)
    
    # Ploting Principal components    
    else:
        pl.plot(pc1,pc2,'k.')
    if title is not None:
        pl.title(title)
    pl.xlabel('PC1')
    pl.xlabel('PC2')
    pl.tight_layout()
    pl.savefig('test',format='pdf')
    pl.clf()
    pl.imshow(kinship_mat, cmap='hot', interpolation='nearest')

    # Ploting the variance explained by each PC 
    tot = sum(evals)
    var_exp = [(i / tot)*100 for i in sorted(evals, reverse=True)]
    np.savetxt('test.out', var_exp, delimiter=',') 
    cum_var_exp = np.cumsum(var_exp)
    print cum_var_exp

    # Ploting the cumulative variance explained
    with pl.style.context('seaborn-whitegrid'):
        pl.figure(figsize=(6, 6))
        pl.bar(range(198), var_exp, alpha= 1, align='center', label='individual explained variance')
        pl.step(range(198), cum_var_exp, where='mid', label='cumulative explained variance')
        pl.ylabel('Explained variance ratio')
        pl.xlabel('Principal components')
        pl.legend(loc='best')
        #plt.tight_layout()
        pl.savefig('variance_explained_PC1234')
        pl.show()

    #pylab.savefig(figure_dir+'/'+k_figure_fn)
plot_dirty_PCA()