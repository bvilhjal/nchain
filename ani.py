import scipy as sp
import h5py 
import pandas as pd
import numpy as np
import pylab as pl
import time
#import matplotlib.pyplot as plt
#import seaborn as sns

def parse_pop_map(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from itertools import izip   
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, origin, country in izip(t['Seq ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'genospecies':origin, 'country':country}
    return pop_map


def compare_equal(array1, array2):
	
	N = len(array1)
	count = 0
	for i in xrange(N):
		if array1[i]==array2[i]:
			count += 1
	ani = count/float(N)

	return ani

def average_nucleotide_identity(snps_hdf5_file = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5', 
                                 seq_file = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/snps.hdf5',
                                 fig_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/figures',
                                 geno_species='gsA', bin_size=0.2,
                                 out_file = '/project/NChain/faststorage/rhizobium/ld/new_snps.hdf5',
                                 gt_hdf5_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/snps.hdf5'):
	pop = parse_pop_map()
	print pop
	pop_map = pop.keys()
	ct_array = pop.values()

	h5f = h5py.File(gt_hdf5_file)
	ag = h5f['alignments']
	gene_big_groups = sorted(ag.keys())
	gene_groups = list()

	# Names
	strains_names = sorted(pop_map, key=lambda x: pop[x]['genospecies'])
	print 'The type of the strain is', type(strains_names)

	strains = []
	for i in ag['1000']['strains']:
		strains.append(i[:4]) 

	print strains
	#map(lambda x: str.replace(x, "3859", "3293"), words)
	print strains_names
	strains_names.remove('3859')
	strains_names.remove('3260')

	print strains_names

	# Making the matrix to update the ani 
	ani_matrix = np.zeros((198,198))
	ani_matrix = pd.DataFrame(ani_matrix, index = strains_names, columns = strains_names)

	# Taking just the core genes
	for gene in gene_big_groups:
		if len(ag[gene]['strains']) == 198:
			gene_groups.append(gene)
	counting = 0

	for gg1 in gene_groups:
		start_time = time.time()
		if counting < 20:
			print 'Working on gene group: %s'%gg1
			print counting
			g1 = ag[gg1]
			for ind1 in zip(g1['nsequences'], g1['strains']): # sequences are in the first entry
				for ind2 in zip(g1['nsequences'], g1['strains']):
					ani_matrix[ind1[1][:4]][ind2[1][:4]] += compare_equal(ind1[0], ind2[0])
		else:
			break
		print "The time used was %6f" % (time.time()-start_time)
		counting += 1	
	ani_matrix = ani_matrix/counting
	
	return ani_matrix
				  
ani = average_nucleotide_identity()
print ani
ani.to_csv('ani_sorted_by_genospecies.csv', header = True)

# A heatmap of the matrix
pl.pcolor(ani)
pl.colorbar()
pl.xlim([0,ani.shape[1]])
pl.ylim([0,ani.shape[0]])
pl.title('Average nucleotide identity - 198 strains')
pl.savefig('heat_map_allcoregenes.png')
pl.show()
