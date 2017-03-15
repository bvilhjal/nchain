# Making LD based on presence absence of 
import itertools 
import numpy as np
from scipy.stats.stats import pearsonr
from collections import Counter
import pandas as pd
import collections

def findInstances(list1, list2):
    """For each item in list1, return a list of offsets to its occurences in list2"""
    for i in list1:
        yield [pos for pos,j in enumerate(list2) if i==j]

# Opening the list of genes and its respective strains
def presence_absence_matrix(genes_directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/PCA_gene_level_all_strains/presence_absence_headers.txt',
                     metadata = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):

  '''It goes from each gene and if a strain have gene it is assigned a number (depending on the genospecies origin) otherwise it is zero (absence of the gene).'''

  genes_strains = open(genes_directory, 'r').read().split()

  # Using ordered dictionary to maintain the gene structure
  gene_groups = collections.OrderedDict() 
  for i in range(len(genes_strains)):
    temp = genes_strains[i].split('|')
    if str(temp[0]) not in gene_groups:
      gene_groups[str(temp[0])] = [int(temp[1])] # creating a list
    else:
      gene_groups[(temp[0])].append(int(temp[1]))

  # Test
  strains_cointained_gene = sorted(gene_groups['group10000'])
  print strains_cointained_gene
    
  # Making the matrix coloured by genospecies

  # Open the data frame with the metadata information
  sorted_strain_names = pd.read_csv(metadata, sep='\t')

  strain_names = list(sorted_strain_names['Seq ID'])
  Matrix_counts = np.zeros((200,len(gene_groups)))
  genospecies = list(sorted_strain_names['Genospecies'])
  countries = {'gsA': 1, 'gsB':1, 'gsC':1, 'gsD':1, 'gsE': 1}

  count = 0
  core = 0
  histogram = list()

  for gene in gene_groups.keys():
    strains = gene_groups[gene]

    index_match = list(findInstances(strains, strain_names))
    if len(index_match) == 199:
      core += 1
    histogram.append(len(index_match))
    merged = list(itertools.chain(*index_match))
    for m in merged:      
      Matrix_counts[m,count] = countries[genospecies[m]]
    count += 1

  print 'Number of genes in total', count
  print 'Number of core genes', core
  print 'The shape of the matrix of counts is:', Matrix_counts.shape

  # Saving the presence absence matrix
  ids = strain_names
  Matrix_counts = np.insert(Matrix_counts, 0, ids, axis=1)
  head = ', '.join(gene_groups.keys())
  head = 'strains,' + head
  np.savetxt("presence_absence_matrix_indexes_1s.csv", Matrix_counts, header = str(head), delimiter=",", fmt='%10.1f')
  return histogram

#presence_absence_matrix()

## Once we have the matrix calculated, we can just download it and do the same calculations of LD
def mantel_test(presence_absence = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/presence_absence_matrix_indexes_1s.csv'):
  pres_abs_matrix = pd.read_csv(presence_absence, sep = ',', header = 0, index_col = 0)
  print pres_abs_matrix.head()

  # Looping through columns (genes)
  for gene1 in pres_abs_matrix:
    for gene2 in pres_abs_matrix:
      print pearsonr(pres_abs_matrix[gene1], pres_abs_matrix[gene2])

mantel_test()

