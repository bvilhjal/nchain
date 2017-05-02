# Creating a venn diagram based on genospecies
import scipy as sp
import itertools 
import matplotlib
import pylab as pl
import numpy as np
import json
import collections
from collections import Counter
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles


def findInstances(list1, list2):
    """For each item in list1, return a list of offsets to its occurences in list2"""
    for i in list1:
        yield [pos for pos,j in enumerate(list2) if i==j]

def parse_pop_map(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, origin, country in izip(t['Seq ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'genospecies':origin, 'country':country}
    return pop_map

def make_histograms(dicts, xlab = None, ylab = None, save_name = None, fig_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/figures'):
  X = np.arange(len(dicts))
  pl.bar(X, sorted(dicts.keys()), align='center', width=0.5)
  pl.xticks(X, sorted(dicts.keys()))
  ymax = max(dicts.values())*1.05
  xmax = len(dicts)
  pl.xlabel(xlab)
  pl.ylabel(ylab)
  pl.ylim(0, ymax)
  pl.xlim(-1,xmax)
  pl.grid(True)
  pl.savefig('%s/figure_%s.pdf'%(fig_dir,save_name))
  pl.show()


#gene_groups = np.load('C:/Users/MariaIzabel/Desktop/MASTER/PHD/Syn-Non/gene_groups/gene_groups.npy').item()

# Opening the list of genes and its respective strains
genes_strains = open('C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/PCA_gene_level_all_strains/presence_absence_headers.txt', 'r').read().split()
# Using ordered dictionary to maintain the gene structure
gene_groups = collections.OrderedDict() 
for i in range(len(genes_strains)):
  temp = genes_strains[i].split('|')
  if str(temp[0]) not in gene_groups:
    gene_groups[str(temp[0])] = [int(temp[1])] # creating a list
  else:
    gene_groups[(temp[0])].append(int(temp[1]))
    

# Making the matrix coloured by genospecies

# Open the data frame with the origin information
sorted_strain_names = pd.read_csv('strains_id.txt', sep='\t')
#print sorted_strain_names
strain_names = list(sorted_strain_names['Seq ID'])
Matrix_counts = np.zeros((200,len(gene_groups)))
lands = list(sorted_strain_names['Genospecies'])
countries = {'gsA': 1, 'gsB':2, 'gsC':3, 'gsD':4, 'gsE': 5}

# Filling up first the core genes
key_sorted_by_length = sorted(gene_groups.values(), key=len)
count = 0

for strains in key_sorted_by_length:
  index_match = list(findInstances(strains, strain_names))
  merged = list(itertools.chain(*index_match))
  for m in merged:      
    Matrix_counts[m,count] = countries[lands[m]]
  count += 1


print Matrix_counts
print 'The shape of the matrix of counts is:', Matrix_counts.shape

# Saving the matrix
np.savetxt("test_presence_absence_matrix.csv", Matrix_counts, delimiter=",")

#Matrix_counts = sorted(Matrix_counts, key=lambda col: sum(col))
#separate core and accessory genes:
plt.pcolor(Matrix_counts[0:100,0:100])
plt.ylim([0,199])
plt.title('Core and accessory genes of 200 strains')
plt.xlim([0,len(gene_groups)])
plt.xlabel('Genes')
plt.ylabel('Strains')
#plt.savefig('presence_absence_genospecies', format = 'pdf')
plt.show()


pop = parse_pop_map()

### Venn diagrams 

def define_combinations(genospecies):
	n = len(genospecies)
	d = {gs:0 for gs in genospecies}

	for k in xrange(2,n+1):
		combinations = itertools.combinations(genospecies, k)
		for comb in combinations:
			d["".join(comb)] = 0
	return d

genospecies = ['gsA', 'gsB', 'gsC', 'gsD', 'gsE']

venn_dictionary = define_combinations(genospecies)
total_geno = {'gsA': 0, 'gsB':0, 'gsC': 0,'gsD': 0, 'gsE': 0}
hist_per_geno = {'gsA': {}, 'gsB':{}, 'gsC': {},'gsD': {}, 'gsE': {}}

for gene in gene_groups:
  strains = gene_groups[gene]
  gs_list = sp.array([pop[str(strain)]['genospecies'] for strain in strains])
  counts_dict = Counter(gs_list)

  # How many genes in total are present in each genospecies?
  gs_list = sp.unique(gs_list)
  for i in gs_list:
    total_geno[i] += 1

  # How many genes contain a given amount of strains from a certain genospecies?:
  for geno, counts in counts_dict.items():
    tup = (geno, counts)
    if tup[1] not in hist_per_geno[tup[0]]:
      hist_per_geno[tup[0]][tup[1]] = 1
    else:
      hist_per_geno[tup[0]][tup[1]] += 1

  # How many genes are shared by each possible combination of genospecies?:    
  gs_list_joined = "".join(sorted(sp.unique(gs_list)))
  venn_dictionary[gs_list_joined] += 1


print '\nAll the possible combinations is:', venn_dictionary
print '\nTotal amount of genes present in each genospecies', total_geno
print '\nDistribution of genes per genospecies', hist_per_geno
json.dump(venn_dictionary, open("text.txt",'w'))


'''Figures'''

### Creating histograms
'''Total number of strains present in each genopecies'''
strains_geno = collections.OrderedDict() 
strains_geno = {'gsC': 116, 'gsB': 34, 'gsA': 33, 'gsE': 11, 'gsD': 5}

'''Total number of genes present in each genospecies'''

make_histograms(total_geno, xlab = 'Genospecies', ylab = 'Genes', save_name = 'genes_genospecies.pdf')

make_histograms(strains_geno, xlab = 'Genospecies', ylab = 'Strains', save_name = 'strains_genospecies.pdf')

for i in strains_geno.keys():
  make_histograms(hist_per_geno[i], xlab= 'Strains', ylab = 'Genes', save_name = i)

#### Creating the venn diagram

# # # Subset sizes
s = (
     1082,   # Abc 
     898,    # aBc 
     332,    # ABc  
     3567,   # abC
     1141,   # AbC
     1290,   # DK/FR
     5904,   # ABC 
  )

v = venn3(subsets=s, set_labels=('Genospecies A', 'Genospecies B', 'Genospecies C'))

# # # Subset labels
v.get_label_by_id('100').set_text('1082')
v.get_label_by_id('010').set_text('898')
v.get_label_by_id('110').set_text('332')
v.get_label_by_id('001').set_text('3567')
v.get_label_by_id('101').set_text('1141')
v.get_label_by_id('011').set_text('1290')
v.get_label_by_id('111').set_text('5904')

# # # Subset colors
v.get_patch_by_id('100').set_color('purple') 	
v.get_patch_by_id('010').set_color('green')
v.get_patch_by_id('110').set_color('orange')

# # # Subset alphas
v.get_patch_by_id('101').set_alpha(0.4)
v.get_patch_by_id('011').set_alpha(1.0)
v.get_patch_by_id('111').set_alpha(0.7)

# # # Border styles
c = venn3_circles(subsets=s, linestyle='solid')
c[0].set_ls('dotted')  # Line style
c[1].set_ls('dashed')
c[2].set_lw(1.0)       # Line width


plt.savefig('venn_diagram_abc.pdf')
plt.clf()


# # ###### Venn Diagram for Genospecies D and E
s = (
   3*83,  # Ab
   3*463,  # aB
   5904,  # AB
	)

v = venn2(subsets=s, set_labels=('Genospecies D', 'Genospecies E'))

# # # Subset labels
v.get_label_by_id('10').set_text('83')
v.get_label_by_id('01').set_text('463')
v.get_label_by_id('11').set_text('5904')

# # # Subset colors
v.get_patch_by_id('10').set_color('c')
v.get_patch_by_id('01').set_color('#993333')
v.get_patch_by_id('11').set_color('blue')

# # # Subset alphas
v.get_patch_by_id('10').set_alpha(0.4)
v.get_patch_by_id('01').set_alpha(1.0)
v.get_patch_by_id('11').set_alpha(0.7)

# # # Border styles
c = venn2_circles(subsets=s, linestyle='solid')
c[0].set_ls('dashed')  # Line style
c[0].set_lw(2.0)       # Line width


plt.savefig('venn_diagram_de.pdf')
plt.clf()
