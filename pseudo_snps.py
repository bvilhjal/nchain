"""
Calculate intergenic LD.
Goal: Identify genes that are more correlated with each other than we expect.
First idea: Correlate gene-specific GRM with each other, accounting for overall covariance (population structure).
6. Identify interesting genes.
7. Write a Nature paper. 
"""

# Cloned the directory to a mac computer
import glob
import numpy as np
import scipy as sp
from scipy import linalg
import h5py
import random
import pandas as pd
import pylab as pl
import time
import simple_intergenic_ld as mantel_test
from sys import argv

# Setting the directory of the data depending on what computer I am  

if argv[1] != 'mac':
    snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'
    out_dir='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps_test/'
    figure_dir='C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Figures/'
else:
    snps_file='/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5'
    out_dir='/Users/PM/Desktop/PHD_incomplete/Methods/Intergenic_LD/corrected_snps_test/'
    figure_dir='/Users/PM/Desktop/PHD_incomplete/nchain/Figures/'    


def calc_ld_table(snps, min_r2 = 0, verbose=True, normalize = True):
    """Calculated LD between all SNPs using r^2 for each gene individually, this function retains snps with values above a given threshold
        Inputed genetic matrix is already normalized by columns/snps"""
    
    if verbose:
        print 'Calculating LD table'

    t0 = time.time()

    # Transposing the snp matrix just to look like Bjarni's function (individuals are in the columns and snps are in the rows)
    snps = snps.T

    # Normalize SNPs (perhaps not necessary, but cheap)
    if normalize:
        snps = (snps - sp.mean(snps, 0)) / sp.std(snps, 0)

    num_snps, num_ind = snps.shape

    # Maximum distance is the distance within a gene
    max_ld_dist = num_snps

    ld_table = {}
    for i in range(num_snps):
        # Initializing the dictionary entries
        ld_table[i] = {}

    a = min(max_ld_dist, num_snps)
    num_pairs = (a * (num_snps - 1)) - a * (a + 1) * 0.5

    if verbose:
        print 'Correlation between %d pairs will be tested' % num_pairs

    num_stored = 0
    for i in range(0, num_snps - 1):
        start_i = i + 1
        end_i = min(start_i + max_ld_dist, num_snps)
        
        ld_vec = sp.dot(snps[i], sp.transpose(snps[start_i:end_i])) / float(num_ind)
        ld_vec = np.array(ld_vec).flatten()

        for k in range(start_i, end_i):
            ld_vec_i = k - start_i
            if ld_vec[ld_vec_i] ** 2 > min_r2:
                ld_table[i][k] = ld_vec[ld_vec_i]
                ld_table[k][i] = ld_vec[ld_vec_i]
                num_stored += 1

    if verbose:
        print 'Stored %d (%0.4f%%) correlations that made the cut r^2 >%0.3f' % (num_stored, 100 * (num_stored / float(num_pairs)), min_r2)

    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60)
    del snps
    return ld_table

def ld_pruning(ld_table, max_ld = 0.99, verbose = True):
    """
    Prunes SNPs in LD, in random order.
    """
    if verbose:
        print 'Prune SNPs in LD, in random order'
    t0 = time.time()
    indices_to_keep = []
    num_snps = len(ld_table)
    indices = sp.random.permutation(num_snps)
    remaining_indices = set(indices)
    for i in indices:
        if len(remaining_indices) == 0:
            break
        elif not (i in remaining_indices):
            continue
        else:
            indices_to_keep.append(i)
            for j in ld_table[i]:
                if ld_table[i][j] > max_ld and j in remaining_indices:
                    remaining_indices.remove(j)

    filter_vector = sp.zeros(num_snps, dtype = 'bool')
    filter_vector[indices_to_keep] = 1
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to LD-prune' % (t / 60, t % 60)
    return filter_vector

def normalize(matrix, direction = 0):
    """Normalizes a matrix default (columns) to mean 0 and var 1."""
    mean = np.mean(matrix, axis= direction)
    std = np.std(matrix, axis= direction)

    matrix = (matrix - mean) / std
    return matrix


def mean_adj(matrix, direction = 1):
    """Normalizes a matrix default (columns) to mean 0 and var 1."""
    mean = np.mean(matrix, axis= direction)
    if direction == 1:
        mean = mean[:, None]

    matrix = (matrix - mean)
    return matrix


def pseudo_snps(snps_file= snps_file,
                 out_dir= out_dir,
                 figure_dir= figure_dir,
                 fig_id='all',
                 min_maf=0,
                 n_snps_delete=10,
                 max_strain_num=200,
                 fig_name='test.png',
                 debug_filter=0.15,
                 write_files = False,
                 min_snps = 5,
                 slicing = True,
                 max_ld = 0.95):

    """
    Take the genes concatenate their snps, calculate GRM, mean adjust the individuals, decomposition, calculate pseudo snps.
    """
    snp_matrices = []
    matrix_lengths = []
    matrix_file_paths = []
    strain_list_masks = []
    snps_to_remove = []

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

    trues = 0
    for i, gg in enumerate(gene_groups):
        if i % 100 == 0:
            print 'Working on gene nr. %d' % i
        if debug_filter >= random.random():
            data_g = h5f[gg]
            strains = data_g['strains'][...]
            if len(strains) < max_strain_num:
                strain_mask = strain_index.get_indexer(strains)

                # Snps M x N
                snps = data_g['norm_snps'][...]

                # Minor allele frequence filtering
                freqs = data_g['freqs'][...]
                mafs = sp.minimum(freqs, 1 - freqs)
                maf_mask = mafs > min_maf
                snps_maf = snps[maf_mask]
                trues += np.sum(maf_mask)

                # Strains in rows and snps in columns:
                snps_maf = snps_maf.T

                if snps_maf.shape[1] >= min_snps:

                    # Changing the precision of the array:
                    snps_maf = np.float64(snps_maf)

                    # LD pruning: Calculate association between samples without over-weighting the contribution of groups of correlated SNPs.
                    ld_table_gene = calc_ld_table(snps_maf, verbose = False)
                    ld_filter_gene = ld_pruning(ld_table_gene, verbose = False, max_ld = max_ld)
                    snps_maf = snps_maf[:,ld_filter_gene]

                    # The SNP matrices are assumed to be sorted by strain. Create a NxM matrix (N = # strains, M = # SNPs) with the
                    # correct rows filled in by the data from the SNP file.
                    full_matrix = np.zeros((198, snps_maf.shape[1]))
                    full_matrix[strain_mask, :] = snps_maf

                    snp_matrices.append(full_matrix)  # the matrix completed with NaNs
                    strain_list_masks.append(strain_mask)  # The strains caring that gene
                    matrix_lengths.append(full_matrix.shape[1])  # The length of the gene
                    matrix_file_paths.append(gg)  # The name of the gene

    snp_boundaries = np.cumsum(matrix_lengths).tolist()

    full_genotype_matrix = np.hstack(snp_matrices)

    print '%d SNPs out of 717457 after filtering for MAF > %f' % (trues, min_maf)
    print 'From the total %d SNPs, %d passed LD pruning' % (trues, full_genotype_matrix.shape[1])

    print 'The full genotype matrix has shape %f' % full_genotype_matrix.shape[1]

    print 'Input the matrix if still there is a nan...'
    full_genotype_matrix = normalize(full_genotype_matrix, direction = 0) 

    print("The variance by column is:...")
    print np.var(full_genotype_matrix, axis = 0)

    print("The mean by columns is:...")
    print np.mean(full_genotype_matrix, axis = 0)

    # 2. Calculate A, the cholesky decomp of the inverse of the GRM.
    print("Finding inverse and sqrt of covariance matrix...")

    N, M = full_genotype_matrix.shape
    not_solved = True

    snp_indices = range(M)

    t0 = time.time()
    while not_solved:

        # Randomly shuffling the SNP columns
        snp_indices_temp = random.sample(snp_indices, len(snp_indices))

        # Making a temporary version of the full matrix
        full_genotype_matrix_temp = full_genotype_matrix[:, snp_indices_temp]

        # Mean adjust the individual lines after shuffling the columns
        print('Normalizing matrix by individuals...')
        full_genotype_matrix_temp = mean_adj(full_genotype_matrix_temp)    

        print('After mean adjust the individuals:')
        print np.mean(full_genotype_matrix_temp, axis = 1)

        # 2. Calculate genome-wide covariance matrix (Kinship)
        print("Calculating genotype matrix covariance...")
        cov = np.dot(full_genotype_matrix_temp, full_genotype_matrix_temp.T) / full_genotype_matrix_temp.shape[1]
        
        try:
            sqrt = linalg.cholesky(cov)
            inv_cov_sqrt = linalg.inv(sqrt.T)
        except:
            continue
        not_solved = False

    t1 = time.time()
    t = (t1 - t0)

    print 'It took %d minutes and %0.2f seconds to solve the inverse of the covariance  matrix' % (t/60, t % 60)
    print inv_cov_sqrt

    # 3. Calculate the pseudo-SNPs (x*A)
    print("Calculating pseudo SNPS...")
    pseudo_snps = np.dot(inv_cov_sqrt, full_genotype_matrix)

    print("The variance by column is:...")
    print np.var(full_genotype_matrix, axis = 0)

    print("The mean by columns is:...")
    print np.mean(full_genotype_matrix, axis = 0)

    print("Dimension")
    print full_genotype_matrix.shape

    del full_genotype_matrix_temp
    del full_genotype_matrix
    
    pl.matshow(cov)
    #pl.title('Kinship - 198 strains - all good genes')
    pl.colorbar()
    pl.savefig(figure_dir + fig_name + 'heat_map_allgenes.pdf')
    pl.show()
  
    identity = np.cov(pseudo_snps)
    pl.matshow(identity)
    #pl.title('After structure correction')
    pl.colorbar()
    pl.savefig(figure_dir + fig_name + 'covariance_pseudo_snps.pdf')
    pl.show()

    print("Creating files corrected for Population Structure...")

    if write_files == True:
    # Extract the original genes from the large pseudo SNP matrix.
        for i, (start, end) in enumerate(zip([0] + snp_boundaries, snp_boundaries)):
            strains_list_mask = strain_list_masks[i]

            if slicing:
                snps = pseudo_snps[strains_list_mask, start:end]
            else:
                snps = pseudo_snps[:,start:end]
            
            strains = strains_list_mask

            file_name = 'group'+matrix_file_paths[i] # the name of the gene

            np.savez_compressed("{}/{}".format(out_dir, file_name), matrix=snps, strains=strains) # structure of the file

#pseudo_snps(min_maf=0.05, fig_name='maf_05_test', debug_filter=1, write_files = True, slicing = True, max_ld = 1)
mantel_test.mantel_corrected_nod_genes(fig_name = 'ld_pruning_1_maf005_with_slicing.pdf', slicing = True)