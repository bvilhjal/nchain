"""
Calculate intergenic LD.

Goal: Identify genes that are more correlated with each other than we expect.

First idea: Correlate gene-specific GRM with each other, accounting for overall GRM (population structure).


6. Identify interesting genes.
7. Write a Nature paper. 


"""
import glob
import bottleneck as bn
import numpy as np
import Mantel
import scipy as sp
from scipy import linalg
import h5py
import random
import pandas as pd
import pylab as pl
import time
#conda install -c https://conda.anaconda.org/biocore scikit-bio

def calc_ld_table(snps, threshold = 0.2, verbose = True, normalize = True):
    """Calculated LD between all SNPs using r^2, this function retains snps with values above a given threshold
        Inputed genetic matrix is already normalized by columns/snps"""
    # Transposing the snp matrix just to look like Bjarni's function
    snps = snps.T 
    num_snps, num_ind = snps.shape

    # Initializing the matrix
    ld_table = np.zeros((num_snps, num_snps))

    # Looping through each comparison, upper triangle of the table
    for i in range(len(ld_table) - 1):
        for j in range(i + 1, len(ld_table)):
            print i



def normalize(matrix, direction = 1):
    """Normalizes a matrix default (columns) to mean 0 and std 1."""
    mean = np.mean(matrix, axis= direction)
    std = np.std(matrix, axis= direction)
    if direction == 1:
        mean = mean[:, None]
        std = std[:, None]

    matrix = (matrix - mean) / std
    return np.nan_to_num(matrix)
            
def replace_column_nans_by_mean(matrix):
    # Set the value of gaps/dashes in each column to be the average of the other values in the column.
    nan_indices = np.where(np.isnan(matrix))
    # Note: bn.nanmean() instead of np.nanmean() because it is a lot(!) faster.
    column_nanmeans = bn.nanmean(matrix, axis=0)

    # For each column, assign the NaNs in that column the column's mean.
    # See http://stackoverflow.com/a/18689440 for an explanation of the following line.
    matrix[nan_indices] = np.take(column_nanmeans, nan_indices[1])

def pseudo_snps(snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5',
                 out_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/',
                 figure_dir='/project/NChain/faststorage/rhizobium/ld/figures',
                 fig_id='all',
                 min_maf=0,
                 max_strain_num=200):
    
    """
    Take the genes concatenate their snps, calculate GRM, decomposition, calculate pseudo snps.
    """
    maf_list_mask = []
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
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        #print strains
        if len(strains) < max_strain_num:
            strain_mask = strain_index.get_indexer(strains)
            snps = data_g['norm_snps'][...]

            freqs = data_g['freqs'][...]
            mafs = sp.minimum(freqs, 1 - freqs)
            maf_mask = mafs >= min_maf
            maf_list_mask.append(maf_mask)
            snps = snps[maf_mask]
            trues += np.sum(maf_mask)

            # Strains in rows and snps in collumns
            snps = snps.T

            # Changing the precision of the array:
            snps = np.float64(snps)

            # The SNP matrices are assumed to be sorted by strain. Create a NxM matrix (N = # strains, M = # SNPs) with the
            # correct rows filled in by the data from the SNP file.
            full_matrix = np.empty((198, snps.shape[1]))
            full_matrix[:] = np.NAN
            full_matrix[strain_mask, :] = snps[:,]

            snp_matrices.append(full_matrix)
    
            strain_list_masks.append(strain_mask)
            matrix_lengths.append(full_matrix.shape[1])
            matrix_file_paths.append(gg) # The name of the gene

    print '%d SNPs out of 717457 after filtering for MAF > %f' % (trues, min_maf)

    snp_boundaries = np.cumsum(matrix_lengths).tolist()
    print("Normalizing genotype matrix...")
    full_genotype_matrix = np.hstack(snp_matrices)

    print('Input the matrix...')
    replace_column_nans_by_mean(full_genotype_matrix)
    #full_genotype_matrix = normalize(full_genotype_matrix, direction=0)
    
    # 2. Calculate A, the cholesky decomp of the inverse of the GRM.
    print("Finding inverse and sqrt of covariance matrix...")

    # Genetically identical individuals results in singular (i.e. non-invertible). This can happen for a subset of the
    # data but should not happen for the entire dataset. Use the pseudo inverse instead. When a matrix is invertible,
    # its pseudo inverse is its inverse anyway. Similarly, pseudo inverses might not be positive definite, meaning that
    # we can't use Cholesky decomposition. If that happens, we can use SciPy's linalg.sqrtm() method (I don't actually
    # know if that is equivalent). Anyway, use linalg.sqrtm(linalg.pinv(cov)) when dealing with small sample sizes.
    
    N, M = full_genotype_matrix.shape
    not_solved = True
    
    shuffle_indices = range(M)
    snp_indices = range(M)

    t0 = time.time()
    while not_solved:
        
        # Random shuffling of the column (SNP) matrix
        #random.shuffle(shuffle_indices)
        #full_genotype_matrix = full_genotype_matrix[:, shuffle_indices]

        # Deleting randomly 10 SNP columns (-10)
        snp_indices_temp = random.sample(snp_indices, len(snp_indices) - 10)

        # Making a temporary version of the full matrix
        full_genotype_matrix_temp = full_genotype_matrix[:, snp_indices_temp]

        # Normalize the individual lines after removing some columns
        print('Normalizing matrix by individuals...')
        full_genotype_matrix_temp = normalize(full_genotype_matrix_temp, direction=1)    

        # Calculate genome-wide GRM/cov (X*X'/M).
        print("Calculating genotype matrix covariance...")
        #cov = np.cov(full_genotype_matrix)

        cov = np.dot(full_genotype_matrix_temp, full_genotype_matrix_temp.T)/full_genotype_matrix_temp.shape[1]

        try:
            inv_cov_sqrt = linalg.cholesky(linalg.inv(cov))
        except:
            continue
        not_solved = False

    t1 = time.time()
    t = (t1 - t0)

    print 'It took %d minutes and %0.2f seconds to solve the inverse of the covariance  matrix' % (t/60, t % 60)
    print inv_cov_sqrt

    # 3. Calculate the pseudo-SNPs (x*A)
    print("Calculating pseudo SNPS...")
    pseudo_snps = np.column_stack(np.dot(inv_cov_sqrt, col) for col in full_genotype_matrix_temp.T)
    del full_genotype_matrix_temp
    del full_genotype_matrix

    pl.matshow(cov)
    pl.title('Kinship - 198 strains - all good genes')
    pl.savefig('heat_map_allgenes.png')

    pl.show()
    pl.matshow(np.cov(pseudo_snps))
    pl.title('After structure correction')
    pl.show()
    pl.savefig('heat_map_structured_genes.png')

    print("Creating files corrected for Population Structure...")

    # Extract the original genes from the large pseudo SNP matrix.
    for i, (start, end) in enumerate(zip([0] + snp_boundaries, snp_boundaries)):
        strains_list_mask = strain_list_masks[i]
        snps = pseudo_snps[strains_list_mask, start:end]
        strains = strains_list_mask
        maf = maf_list_mask[i]

        file_name = 'group'+matrix_file_paths[i] # the name of the gene
        
        np.savez_compressed("{}/{}".format(out_dir, file_name), matrix=snps, strains=strains, maf = maf) # structure of the file
    
print pseudo_snps()


def intergenic_ld(in_glob = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/*.npz'):
    distance_to_ld = {}

    #This is the gene SNPs matrix
    genes = []
    for f in glob.glob("*.npz"):
        with np.load(f) as data:
            # Creating a tuple
            genes.append((data["matrix"], data["strains"])) 
    
    r_scores = []
    p_values = []
    z_scores = []
    count = 0
    for gene1 in genes:
        print gene1
        for gene2 in genes:
            if count <= 1000:
                print 'working on %d', (count)
            # 5. Perform Mantel test between all genes
                if gene1[0].shape[0] == gene2[0].shape[0]:
                    g1 = np.dot(gene1[0], gene1[0].T)/gene1[0].shape[1]
                
                    # Tranforming the diagonal to zero
                    g1 = g1 - np.diag(np.diag(g1))

                    g2 = np.dot(gene2[0], gene2[0].T)/gene2[0].shape[1]
                    g2 = g2 - np.diag(np.diag(g2))
                    r, p, z = Mantel.mantel_test(g1, g2, perms=10, method='pearson', tail='two-tail')
                    count += 1
                    r_scores.append(r)
                    p_values.append(p)
                    z_scores.append(z)

    LD_stats = pd.DataFrame(
    {'r_scores': r_scores,
    'p_values': p_values,
    'z_scores': z_scores
    })
    LD_stats.to_csv('intergenic_LD_stats.csv', header = True)
    return LD_stats

#print intergenic_ld()

