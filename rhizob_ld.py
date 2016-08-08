"""
Code for analysing LD of rhizobium
"""

import scipy as sp
import h5py 
import os

def gen_genotype_hdf5file(out_hdf5_file ='/faststorage/project/NChain/bjarni/snps.hdf5', 
                          directory='/faststorage/project/NChain/rhizobium/ld/snps/'):
    snps_files = os.listdir(directory)
    h5f = h5py.File(out_hdf5_file)
    g = h5f.create_group('gene_groups')
    for i, snps_f in enumerate(snps_files):
        if i%100==0:
            print i
        l = snps_f.split('.')
        group_num = int(l[0][5:])
        snps_data = sp.load(directory+'/'+snps_f)
        snps_g = g.create_group('%d'%group_num)
        snps_g.create_dataset('snps', data=snps_data['matrix'], compression='lzf')
        snps_g.create_dataset('alignment_length',data=snps_data['alignment_length'])
        snps_g.create_dataset('minor_frequencies',data=snps_data['minor_frequencies'])
        snps_g.create_dataset('strains',data=snps_data['strains'].tolist())
        snps_g.create_dataset('positions',data=snps_data['positions'])
        h5f.flush()
    h5f.close()
    