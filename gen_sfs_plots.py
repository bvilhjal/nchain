###### Re-do the plots now with more genospecies:
"""
Code for analysing LD of rhizobium
"""
import numpy as np
import scipy as sp
import h5py 
import os
import csv
import matplotlib
matplotlib.use('Agg')
import pylab
from matplotlib import pyplot as plt
# pylab.rcParams['legend.numpoints'] = 1

import collections
import pandas as pd
import cPickle
import gzip
from sys import argv

def parse_pop_map(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, origin, country in izip(t['Seq ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'genospecies':origin, 'country':country}
    return pop_map

def gen_sfs_plots(snps_hdf5_file = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/called_snps.hdf5', 
                 fig_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/', filter_pop=None):
    
    ### Here I will do the SFS for each genospecies based on the rhizobium xls file
    pop = parse_pop_map()
    pop_map = pop.keys()
    ct_array = pop.values()
    from itertools import izip
    h5f = h5py.File(snps_hdf5_file)
    gene_groups = sorted(h5f.keys())
    
    syn_mafs = []
    nonsyn_mafs = []
    all_mafs = []
    sfs_dict = {}
    for i, gg in enumerate(gene_groups):
        if i%100==0:
            print '%d: Gene %s'%(i,gg)  
        g = h5f[gg]
        if g['codon_snps'].size>1:
            #print g['codon_snps'].shape
            
            if filter_pop is not None:
                strains = g['strains']
                indiv_filter = sp.zeros((len(strains)),dtype='bool8')
                for s_i, s in enumerate(strains):
                    if pop[s]['genospecies']== filter_pop:
                        indiv_filter[s_i]=True
                    codon_snps = g['codon_snps'][...]
                    codon_snps = codon_snps[:,indiv_filter] # reducing the collumns based on the genospecies
                    t_codon_snps = sp.transpose(codon_snps)
                    freqs = sp.mean(t_codon_snps,0)
                    # rows are snps collumns are individuals  
                    #counts = np.sum(codon_snps, axis = 0)
                    #print counts
                    #for c in counts:
                    #    if c in sfs_dict:
                    #        sfs_dict[c] += 1
                    #    else:
                    #        sfs_dict[c] = 1
                #with open('dict.csv', 'wb') as csv_file:
                #    writer = csv.writer(csv_file)
                #    for key, value in sfs_dict.items():
                #        writer.writerow([key, value])
                
            else:
                codon_snps = g['codon_snps'][...]
                t_codon_snps = sp.transpose(codon_snps)
                freqs = sp.mean(t_codon_snps,0) # number of minor allele
            mafs = sp.minimum(freqs,1-freqs)
            is_synonimous_snp = g['is_synonimous_snp'][...]
            syn_mafs.extend(mafs[is_synonimous_snp])
            nonsyn_mafs.extend(mafs[sp.negative(is_synonimous_snp)])
            all_mafs.extend(mafs)
             
    if filter_pop is not None:
        output_file = "%s.csv" %(str(argv[1]))
        np.savetxt(output_file, all_mafs, delimiter=',')   # X is an array
        output_file = "%ssyn_mafs.csv" %(str(argv[1]))
        np.savetxt(output_file, syn_mafs, delimiter=',')
        output_file = "%snon_syn_mafs.csv" %(str(argv[1])) 
        np.savetxt(output_file, nonsyn_mafs, delimiter=',') 
       # pylab.clf()
       # pylab.hist(all_mafs, bins=50)
       # pylab.title('SFS (all binary codon SNPs)')
       # pylab.savefig('%s/sfs_all_%s.png'%(fig_dir,filter_pop))
    
        #pylab.clf()
        #pylab.hist(nonsyn_mafs, bins=50)
        #pylab.title('SFS (non-synonimous SNPs)')
        #pylab.savefig('%s/sfs_non_syn_%s.png'%(fig_dir,filter_pop))
    
        #pylab.clf()
        #pylab.hist(syn_mafs, bins=50)
        #pylab.title('SFS (synonimous SNPs)')
        #pylab.savefig('%s/sfs_syn_%s.png'%(fig_dir,filter_pop))
        
    else:
        output_file = "total_2.csv"
        np.savetxt(output_file, all_mafs, delimiter=',')
        pylab.clf()
        pylab.hist(all_mafs, bins=50)
        pylab.title('SFS (all binary codon SNPs)')
        pylab.savefig(fig_dir+'/sfs_all.png')
    
        pylab.clf()
        pylab.hist(nonsyn_mafs, bins=50)
        pylab.title('SFS (non-synonimous SNPs)')
        pylab.savefig(fig_dir+'/sfs_non_syn.png')
    
        pylab.clf()
        pylab.hist(syn_mafs, bins=50)
        pylab.title('SFS (synonimous SNPs)')
        pylab.savefig(fig_dir+'/sfs_syn.png')
        
      
gen_sfs_plots(snps_hdf5_file = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/called_snps.hdf5', 
                 fig_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/', filter_pop=argv[1])