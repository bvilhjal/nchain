"""
Code for analysing LD of rhizobium
"""

import scipy as sp
import h5py 
import os

def gen_genotype_hdf5file(out_hdf5_file ='/faststorage/project/NChain/bjarni/snps.hdf5', 
                          snps_directory='/faststorage/project/NChain/rhizobium/ld/snps/',
                          snps_wo_struct_directory='/faststorage/project/NChain/rhizobium/ld/snps_no_structure/'):
    snps_files = os.listdir(snps_directory)

    print "Raw SNPs files"
    h5f = h5py.File(out_hdf5_file)
    g = h5f.create_group('gene_groups')
    for i, snps_f in enumerate(snps_files):
        if i%100==0:
            print i
        l = snps_f.split('.')
        group_num = int(l[0][5:])
        snps_data = sp.load(snps_directory+'/'+snps_f)
        snps_g = g.create_group('%d'%group_num)
        snps_g.create_dataset('snps', data=snps_data['matrix'], compression='lzf')
        snps_g.create_dataset('alignment_length',data=snps_data['alignment_length'])
        snps_g.create_dataset('minor_frequencies',data=snps_data['minor_frequencies'])
        snps_g.create_dataset('strains',data=[a.encode('utf8') for a in snps_data['strains']])
        snps_g.create_dataset('positions',data=snps_data['positions'])
        h5f.flush()

    print "Now SNPs wo structure"
    snps_files = os.listdir(snps_wo_struct_directory)
    for i, snps_f in enumerate(snps_files):
        if i%100==0:
            print i
        l = snps_f.split('.')
        group_num = int(l[0][5:])
        snps_data = sp.load(snps_wo_struct_directory+'/'+snps_f)
        snps_g = g['%d'%group_num]
        snps_g.create_dataset('snps_wo_struct', data=snps_data['matrix'], compression='lzf')        
    h5f.close()
    