"""
Calculate intergenic LD.

Goal: Identify genes that are more correlated with each other than we expect.

First idea: Correlate gene-specific GRM with each other, accounting for overall GRM (population structure).


6. Identify interesting genes.
7. Write a Nature paper. 


"""
import kinship
import scipy as sp
from scipy.sparse import identity
from scipy import linalg
import h5py
import pandas as pd
from skbio.math.stats.distance import mantel

def normalize(matrix, direction = 1):
    """Normalizes a matrix default (columns) to mean 0 and std 1."""
    mean = np.mean(matrix, axis= direction)
    std = np.std(matrix, axis= direction)
    matrix = (matrix - mean) / std
    return np.nan_to_num(matrix)

def pseudo_snps_GRM(matrix):
    """ Calculate the pseudo gene matrix and returns the GRM"""
    standard_matrix = normalize(matrix, direction = 1)
    chol = linalg.cholesky(linalg.pinv(standard_matrix))  
    pseudo_snps = sp.dot(chol, standard_matrix)
    GRM = sp.dot(pseudo_snps.T, pseudo_snps)/pseudo_snps.shape[1]
    return GRM


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
    
    # 1. Calculate genome-wide GRM (X*X'/M).
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
    pseudo_snps = cholesky_decomp_inv_snp_cov*K_snps
    GRM = sp.dot(pseudo_snps.T, pseudo_snps)/pseudo_snps.shape[0] # dont remember 

    # 4. Use pseudo-SNPs to calculate gene-specific GRMs.
    # Here we calculate the gene-GRM

    # 5. Perform Mantel test between all genes.
    # How do we deal with matrices of different shapes?
    # Genes correlating with GRM are respecting genospecies boundaries

    # Genes that do not correlates with GRM are disrespecting genospecies boundaries and may be functionally important for all the genospecies (like symbiotic genes)

    # To simplify, let's first look at genes with the same dimensions
    core_genes = []
    for g in gene_groups:
        if len(gene_groups) == 198:
            core_genes.append(g)

    for g1 in core_genes:
        for g2 in core_genes:
            if g1 != g2:
                GRM_g1 = pseudo_snps_GRM(g1)
                GRM_g2 = pseudo_snps_GRM(g2)

                # Mantel test
                coeff, p_value = mantel(GRM_g1, GRM_g2, method='pearson', permutations=100, alternative='two-sided')





    
    
    
    
    
    
    
