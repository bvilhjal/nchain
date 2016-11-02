"""
Code to calculate the kinship.
"""

import scipy as sp
import h5py
import pandas as pd
import numpy as np 
import matplotlib
import rhizob_ld
matplotlib.use('Agg')
import pylab

    
def get_kinships(snps_file='/project/NChain/faststorage/rhizobium/ld/new_snps.hdf5',
                 plot_figures=False,
                 figure_dir='/project/NChain/faststorage/rhizobium/ld/figures',
                 fig_id='all',
                 min_maf=0.1,
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
    K_codon_snps = sp.zeros((num_strains, num_strains))
    counts_mat_codon_snps = sp.zeros((num_strains, num_strains))
        
    K_nonsyn_snps = sp.zeros((num_strains, num_strains))
    counts_mat_nonsyn_snps = sp.zeros((num_strains, num_strains))
    
    K_syn_snps = sp.zeros((num_strains, num_strains))
    counts_mat_syn_snps = sp.zeros((num_strains, num_strains))

    for i, gg in enumerate(gene_groups):
        if i % 100 == 0:
            print 'Working on gene nr. %d' % i 
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) < max_strain_num:
            strain_mask = strain_index.get_indexer(strains)
            
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
    
            codon_snps = data_g['norm_codon_snps'][...]
            if len(codon_snps) == 0:
                continue
            freqs = data_g['codon_snp_freqs'][...]
            mafs = sp.minimum(freqs, 1 - freqs)
            maf_mask = mafs > min_maf
            codon_snps = codon_snps[maf_mask]
            is_synonimous_snp = data_g['is_synonimous_snp'][...]
            is_synonimous_snp = is_synonimous_snp[maf_mask]
            if len(codon_snps) > 0:
                K_codon_snps_slice = K_codon_snps[strain_mask]
                K_codon_snps_slice[:, strain_mask] += sp.dot(codon_snps.T, codon_snps)
                K_codon_snps[strain_mask] = K_codon_snps_slice
                counts_mat_codon_snps_slice = counts_mat_codon_snps[strain_mask]
                counts_mat_codon_snps_slice[:, strain_mask] += len(codon_snps)
                counts_mat_codon_snps[strain_mask] = counts_mat_codon_snps_slice
        
                if sp.sum(is_synonimous_snp) > 0:
                    syn_snps = codon_snps[is_synonimous_snp]
                    K_syn_snps_slice = K_syn_snps[strain_mask]
                    K_syn_snps_slice[:, strain_mask] += sp.dot(syn_snps.T, syn_snps)
                    K_syn_snps[strain_mask] = K_syn_snps_slice
                    counts_mat_syn_snps_slice = counts_mat_syn_snps[strain_mask]
                    counts_mat_syn_snps_slice[:, strain_mask] += len(syn_snps)
                    counts_mat_syn_snps[strain_mask] = counts_mat_syn_snps_slice
            
                is_nonsynonimous_snp = sp.negative(is_synonimous_snp)
                if sp.sum(is_nonsynonimous_snp) > 0:
                    nonsyn_snps = codon_snps[is_nonsynonimous_snp]                
                    K_nonsyn_snps_slice = K_nonsyn_snps[strain_mask]
                    K_nonsyn_snps_slice[:, strain_mask] += sp.dot(nonsyn_snps.T, nonsyn_snps)
                    K_nonsyn_snps[strain_mask] = K_nonsyn_snps_slice
                    counts_mat_nonsyn_snps_slice = counts_mat_nonsyn_snps[strain_mask]
                    counts_mat_nonsyn_snps_slice[:, strain_mask] += len(nonsyn_snps)
                    counts_mat_nonsyn_snps[strain_mask] = counts_mat_nonsyn_snps_slice

    
    
    K_snps = K_snps / counts_mat_snps  # element-wise division
    K_codon_snps = K_codon_snps / counts_mat_codon_snps  # element-wise division

    K_syn_snps = K_syn_snps / counts_mat_syn_snps  # element-wise division
    K_nonsyn_snps = K_nonsyn_snps / counts_mat_nonsyn_snps  # element-wise division

    if plot_figures:
        plot_dirty_PCA(K_snps, figure_fn='PCA_all_snps_%s.pdf' % fig_id, k_figure_fn='K_all_snps_%s.png' % fig_id,
                       figure_dir=figure_dir, strains=ordered_strains, title='All SNPs')
        plot_dirty_PCA(K_codon_snps, figure_fn='PCA_codon_snps_%s.pdf' % fig_id, k_figure_fn='K_codon_snps_%s.png' % fig_id,
                       figure_dir=figure_dir, strains=ordered_strains, title='Codon SNPs')
        plot_dirty_PCA(K_syn_snps, figure_fn='PCA_syn_snps_%s.pdf' % fig_id, k_figure_fn='K_syn_snps_%s.png' % fig_id,
                       figure_dir=figure_dir, strains=ordered_strains, title='Synonymous SNPs')
        plot_dirty_PCA(K_nonsyn_snps, figure_fn='PCA_nonsyn_snps_%s.pdf' % fig_id, k_figure_fn='K_nonsyn_snps_%s.png' % fig_id,
                       figure_dir=figure_dir, strains=ordered_strains, title='Non-Synonymous SNPs')

    print 'Average number of SNPs: %0.2f.' % sp.mean(counts_mat_snps)
    print 'Average number of codon SNPs: %0.2f.' % sp.mean(counts_mat_snps)
    print 'Average number of codon SNPs: %0.2f.' % sp.mean(counts_mat_snps)
    print 'Average number of codon SNPs: %0.2f.' % sp.mean(counts_mat_snps)

    return {'K_snps':K_snps, 'K_codon_snps':K_codon_snps, 'counts_mat_snps':counts_mat_snps, 'counts_mat_codon_snps':counts_mat_codon_snps,
            'K_syn_snps':K_syn_snps, 'K_nonsyn_snps':K_nonsyn_snps, 'counts_mat_syn_snps':counts_mat_syn_snps, 'counts_mat_nonsyn_snps':counts_mat_nonsyn_snps,
            'strains':ordered_strains}
    

def plot_dirty_PCA(kinship_mat, figure_fn='pca.png', k_figure_fn='kinship_heatmap.png', title=None,
                   figure_dir='/project/NChain/faststorage/rhizobium/ld/figures', strains=None):
    from scipy import linalg
    
    evals, evecs = linalg.eig(kinship_mat)  # PCA via eigen decomp
    evals[evals < 0] = 0
    sort_indices = sp.argsort(evals,)
    ordered_evals = evals[sort_indices]
    print ordered_evals[-10:] / sp.sum(ordered_evals)
    pc1, pc2 = evecs[:, sort_indices[-1]], evecs[:, sort_indices[-2]]
    pylab.clf()
    
    
    if strains is not None:    
        ct_marker_map = {'DK':'*', 'UK':'^', 'F':'o', 'NA': 's'}
        gs_color_map = {'gsA':'m', 'gsB':'g', 'gsC':'r', 'gsE': 'b', 'NA':'c'}
        pop_map = rhizob_ld.parse_pop_map()
        for i, strain in enumerate(strains):
            d = pop_map.get(strain, 'NA')
            if d == 'NA':
                gs = 'NA'
                country = 'NA'
            else:
                gs = d['genospecies']
                country = d['country']
            pylab.scatter(pc1[i], pc2[i], marker=ct_marker_map[country], c=gs_color_map[gs], alpha=0.3, s=100, edgecolor='none')
        for gs in gs_color_map:
            pylab.scatter([], [], color=gs_color_map[gs], marker='s', label=gs, s=100, edgecolor='none')
        for country in ct_marker_map:
            if country != 'NA':
                pylab.scatter([], [], color='k', marker=ct_marker_map[country], label=country, s=100, facecolors='none')

        
        pylab.legend(scatterpoints=1)
        
    else:
        pylab.plot(pc1, pc2, 'k.')
    if title is not None:
        pylab.title(title)
    pylab.xlabel('PC1')
    pylab.xlabel('PC2')
    pylab.tight_layout()
    pylab.savefig(figure_dir + '/' + figure_fn, format='pdf')
    pylab.clf()
    pylab.imshow(kinship_mat, cmap='hot', interpolation='nearest')
    pylab.savefig(figure_dir + '/' + k_figure_fn)

    