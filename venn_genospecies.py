# Creating a venn diagram based on genospecies

import scipy as sp
import itertools 
import matplotlib
import pylab as pl
import numpy as np
import json
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles

def parse_pop_map(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, origin, country in izip(t['Seq ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'genospecies':origin, 'country':country}
    return pop_map

gene_groups = np.load('C:/Users/MariaIzabel/Desktop/MASTER/PHD/Syn-Non/gene_groups/gene_groups.npy').item()

pop = parse_pop_map()

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

for gene in gene_groups:
	strains = gene_groups[gene]
	gs_list = sp.array([pop[strain]['genospecies'] for strain in strains])
	gs_list = "".join(sorted(sp.unique(gs_list)))
	#print gs_list
	venn_dictionary[gs_list] += 1

print 'All the possible combinations is: ', venn_dictionary

json.dump(venn_dictionary, open("text.txt",'w'))

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
