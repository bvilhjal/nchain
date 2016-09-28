###### Re-do the plots now with more genospecies:
"""
Code for analysing LD of rhizobium
"""
import numpy as np
import scipy as sp
import h5py 
import os
import matplotlib
matplotlib.use('Agg')
import pylab
from matplotlib import pyplot as plt
# pylab.rcParams['legend.numpoints'] = 1

import collections
import pandas as pd
import cPickle
import gzip

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
    for i, gg in enumerate(gene_groups):
        if i%100==0:
            print '%d: Gene %s'%(i,gg)  
        g = h5f[gg]
        if g['codon_snp_freqs'].size>1:
            
            if filter_pop is not None:
                strains = g['strains']
                indiv_filter = sp.zeros((len(strains)),dtype='bool8')
                for s_i, s in enumerate(strains):
                    try:
                        if pop[s]['genospecies']== filter_pop:
                            indiv_filter[s_i]=True
                    except:
                        continue
                if sp.sum(indiv_filter)<100:
                    continue
                codon_snps = g['codon_snps'][...]
                codon_snps = codon_snps[:,indiv_filter]
                t_codon_snps = sp.transpose(codon_snps)
                freqs = sp.mean(t_codon_snps,0)
                #print freqs
                
            else:
                freqs = g['codon_snp_freqs'][...]
            mafs = sp.minimum(freqs,1-freqs)
            is_synonimous_snp = g['is_synonimous_snp'][...]
            syn_mafs.extend(mafs[is_synonimous_snp])
            nonsyn_mafs.extend(mafs[sp.negative(is_synonimous_snp)])
            all_mafs.extend(mafs)
    np.savetxt('test_A.csv', all_mafs, delimiter=',')   # X is an array
            
    
    if filter_pop is not None:
        pylab.clf()
        pylab.hist(all_mafs, bins=50)
        pylab.title('SFS (all binary codon SNPs)')
        pylab.savefig('%s/sfs_all_%s.png'%(fig_dir,filter_pop))
    
        #pylab.clf()
        #pylab.hist(nonsyn_mafs, bins=50)
        #pylab.title('SFS (non-synonimous SNPs)')
        #pylab.savefig('%s/sfs_non_syn_%s.png'%(fig_dir,filter_pop))
    
        #pylab.clf()
        #pylab.hist(syn_mafs, bins=50)
        #pylab.title('SFS (synonimous SNPs)')
        #pylab.savefig('%s/sfs_syn_%s.png'%(fig_dir,filter_pop))
        
    else:
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
                 fig_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/', filter_pop='gsA')