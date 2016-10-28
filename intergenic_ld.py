"""
Calculate intergenic LD.

Goal: Identify genes that are more correlated with each other than we expect.

First idea: Correlate gene-specific GRM with each other, accounting for overall GRM (population structure).


4. Use pseudo-SNPs to calculate gene-specific GRMs.
5. Perform Mantel test between all genes.
6. Identify intersting genes.
7. Write a Nature paper. 


"""
import kinship
import scipy as sp
from scipy import linalg
import h5py
import pandas as pd

def get_snps(snps_file='/project/NChain/faststorage/rhizobium/ld/new_snps.hdf5',
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
            
            

if __name__ == '__main__':
    
    # 1. Calcuate genome-wide GRM.
    k_dict = kinship.get_kinships(snps_file='/project/NChain/faststorage/rhizobium/ld/new_snps.hdf5',
                 plot_figures=False,
                 figure_dir='/project/NChain/faststorage/rhizobium/ld/figures',
                 fig_id='all',
                 min_maf=0.1,
                 max_strain_num=200)
    
    K_snps = k_dict['K_snps']
    
    # 2. Calculate A, the cholesky decomp of the inverse of the GRM.
    # Approximation.. may be improved by using a better SNP covariance matrix.
    cholesky_decomp_inv_snp_cov = linalg.cholesky(linalg.pinv(K_snps))  
    
    # 3. Calculate the pseudo-SNPs (x*A)
    
    
    
    
    
    
    
