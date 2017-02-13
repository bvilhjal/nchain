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


def minor_allele_filter(gene_matrix, maf):
    '''Filtering for minor allele frequency, it assumes that the matrix is binary and that is n x m, where columns are markers (m).'''

    freqs = sp.mean(gene_matrix, 0) 
    mafs = sp.minimum(freqs, 1 - freqs)

    maf_filter = mafs > maf
    matrix_mafs = gene_matrix[:, maf_filter]

    # Normalizing the columns:
    norm_matrix = (matrix_mafs - np.mean(matrix_mafs, axis=0)) / np.std(matrix_mafs, axis=0)

    # Mean adjuste by rows:
    # ???

    return(norm_matrix)



def mantel_test(X, Y, perms=10000, method='pearson', tail='two-tail'):
  """
  Takes two distance matrices (either redundant matrices or condensed vectors)
  and performs a Mantel test. The Mantel test is a significance test of the
  correlation between two distance matrices.

  Parameters
  ----------
  X : array_like
      First distance matrix (condensed or redundant).
  Y : array_like
      Second distance matrix (condensed or redundant), where the order of
      elements corresponds to the order of elements in the first matrix.
  perms : int, optional
      The number of permutations to perform (default: 10000). A larger number
      gives more reliable results but takes longer to run. If the actual number
      of possible permutations is smaller, the program will enumerate all
      permutations. Enumeration can be forced by setting this argument to 0.
  method : str, optional
      Type of correlation coefficient to use; either 'pearson' or 'spearman'
      (default: 'pearson').
  tail : str, optional
      Which tail to test in the calculation of the empirical p-value; either
      'upper', 'lower', or 'two-tail' (default: 'two-tail').

  Returns
  -------
  r : float
      Veridical correlation
  p : float
      Empirical p-value
  z : float
      Standard score (z-score)
  """

  # Ensure that X and Y are formatted as Numpy arrays.
  X, Y = np.asarray(X, dtype=float), np.asarray(Y, dtype=float)

  # Check that X and Y are valid distance matrices.
  # if spatial.distance.is_valid_dm(X) == False and spatial.distance.is_valid_y(X) == False:
  #  raise ValueError('X is not a valid condensed or redundant distance matrix')
  # if spatial.distance.is_valid_dm(Y) == False and spatial.distance.is_valid_y(Y) == False:
  #  raise ValueError('Y is not a valid condensed or redundant distance matrix')

  # If X or Y is a redundant distance matrix, reduce it to a condensed distance matrix.
  if len(X.shape) == 2:
    X = spatial.distance.squareform(X, force='tovector', checks=False)
  if len(Y.shape) == 2:
    Y = spatial.distance.squareform(Y, force='tovector', checks=False)

  # Check for size equality.
  if X.shape[0] != Y.shape[0]:
    raise ValueError('X and Y are not of equal size')

  # Check for minimum size.
  if X.shape[0] < 3:
    raise ValueError('X and Y should represent at least 3 objects')

  # If Spearman correlation is requested, convert X and Y to ranks.
  if method == 'spearman':
    X, Y = stats.rankdata(X), stats.rankdata(Y)

  # Check for valid method parameter.
  elif method != 'pearson':
    raise ValueError('The method should be set to "pearson" or "spearman"')

  # Check for valid tail parameter.
  if tail != 'upper' and tail != 'lower' and tail != 'two-tail':
    raise ValueError('The tail should be set to "upper", "lower", or "two-tail"')

  # Now we're ready to start the Mantel test using a number of optimizations:
  #
  # 1. We don't need to recalculate the pairwise distances between the objects
  #    on every permutation. They've already been calculated, so we can use a
  #    simple matrix shuffling technique to avoid recomputing them. This works
  #    like memoization.
  #
  # 2. Rather than compute correlation coefficients, we'll just compute the
  #    covariances. This works because the denominator in the equation for the
  #    correlation coefficient will yield the same result however the objects
  #    are permuted, making it redundant. Removing the denominator leaves us
  #    with the covariance.
  #
  # 3. Rather than permute the Y distances and derive the residuals to calculate
  #    the covariance with the X distances, we'll represent the Y residuals in
  #    the matrix and shuffle those directly.
  #
  # 4. If the number of possible permutations is less than the number of
  #    permutations that were requested, we'll run a deterministic test where
  #    we try all possible permutations rather than sample the permutation
  #    space. This gives a faster, deterministic result.

  # Calculate the X and Y residuals, which will be used to compute the
  # covariance under each permutation.
  X_residuals, Y_residuals = X - X.mean(), Y - Y.mean()

  # Expand the Y residuals to a redundant matrix.
  Y_residuals_as_matrix = spatial.distance.squareform(Y_residuals, force='tomatrix', checks=False)

  # Get the number of objects.
  m = Y_residuals_as_matrix.shape[0]

  # Calculate the number of possible matrix permutations.
  n = np.math.factorial(m)

  # Initialize an empty array to store temporary permutations of Y_residuals.
  Y_residuals_permuted = np.zeros(Y_residuals.shape[0], dtype=float)

  # If the number of requested permutations is greater than the number of
  # possible permutations (m!) or the perms parameter is set to 0, then run a
  # deterministic Mantel test ...
  if perms >= n or perms == 0:

    # Initialize an empty array to store the covariances.
    covariances = np.zeros(n, dtype=float)

    # Enumerate all permutations of row/column orders and iterate over them.
    for i, order in enumerate(permutations(range(m))):

      # Take a permutation of the matrix.
      Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      spatial.distance._distance_wrap.to_vector_from_squareform_wrap(Y_residuals_as_matrix_permuted, Y_residuals_permuted)

      # Compute and store the covariance.
      covariances[i] = (X_residuals * Y_residuals_permuted).sum()

  # ... otherwise run a stochastic Mantel test.
  else:

    # Initialize an empty array to store the covariances.
    covariances = np.zeros(perms, dtype=float)

    # Initialize an array to store the permutation order.
    order = np.arange(m)

    # Store the veridical covariance in 0th position...
    covariances[0] = (X_residuals * Y_residuals).sum()

    # ...and then run the random permutations.
    for i in range(1, perms):

      # Choose a random order in which to permute the rows and columns.
      np.random.shuffle(order)

      # Take a permutation of the matrix.
      Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      spatial.distance._distance_wrap.to_vector_from_squareform_wrap(Y_residuals_as_matrix_permuted, Y_residuals_permuted)

      # Compute and store the covariance.
      covariances[i] = (X_residuals * Y_residuals_permuted).sum()
      
  # Calculate the veridical correlation coefficient from the veridical covariance.
  r = covariances[0] / np.sqrt((X_residuals ** 2).sum() * (Y_residuals ** 2).sum())

  # Calculate the empirical p-value for the upper or lower tail.
  if tail == 'upper':
    p = (covariances >= covariances[0]).sum() / float(covariances.shape[0])
  elif tail == 'lower':
    p = (covariances <= covariances[0]).sum() / float(covariances.shape[0])
  elif tail == 'two-tail':
    p = (abs(covariances) >= abs(covariances[0])).sum() / float(covariances.shape[0])

  # Calculate the standard score.
  z = (covariances[0] - covariances.mean()) / covariances.std()

  return r, p, z


def simple_intergenic_ld_core(max_strain_num=198,
                            maf=0.2,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    all_strains = set()
    
    core_genes = []
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(set(strains)) == max_strain_num:
            core_genes.append(gg)
    
    ordered_strains = sorted(list(all_strains))
    strain_index = pd.Index(ordered_strains)

    count = 0
    r_scores = []
    p_values = []
    z_scores = []
    gene_names = []

    gene_grm_dict = {}

    for i, gg1 in enumerate(core_genes):
        data_g1 = h5f[gg1]
        total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns 
        print total_snps_1.shape
        total_snps_1 = minor_allele_filter(total_snps_1, 0.1)

        ''' 3. Calculate the Kinship matrix for each gene '''
        grm = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
        print grm.shape
        print np.average(np.diag(grm))

        flat_grm1 = grm.flatten()
        gene_grm_dict[gg1] = {'grm':grm , 'flat_grm':flat_grm1 - flat_grm1.mean()}
        
        for j, gg2 in enumerate(core_genes):
            if i > j:
            
                data_g2 = h5f[gg2]  # tuple

                # Take the subset of shared snps of each gene
                total_snps_2 = data_g2['snps'][...].T
                print total_snps_2.shape


                '''Filtering for Minor allele frequency'''
                total_snps_2 = minor_allele_filter(total_snps_2, 0.1)
                
                flat_grm2 = gene_grm_dict[gg2]['flat_grm']
                r = sp.dot(flat_grm1, flat_grm2) / sp.sqrt(sp.dot(flat_grm1, flat_grm1), sp.dot(flat_grm2, flat_grm2))
                print r
                r_scores.append(r)
                gene_names.append(gg1 + gg2)

    LD_stats = pd.DataFrame(
    {'r_scores': r_scores,
    'p_values': p_values,
    'z_scores': z_scores,
    'gene_names': gene_names})

    LD_stats.to_csv('test.csv', header=True)


    return LD_stats

simple_intergenic_ld_core()


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

