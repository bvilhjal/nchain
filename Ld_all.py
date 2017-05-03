# LD of all 
# 1. Find a reference genome to sort the genes by location
# 2. Calculate the LD of all versus all
# 3. Find possible LD blocks, besides the one we know (nod genes)
# 4. Find gene functions
# 1. Find a reference genome to sort the genes by location

import numpy as np
from collections import OrderedDict
import os
import glob 
import pandas as pd
import scipy as sp
from itertools import izip
from scipy.stats.stats import pearsonr

def gene_locations(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/gene_locations_sorted.csv'):
    '''Adding the location of the genes based on a blast of the pacbio genome SM158'''
    
    loc_map = {}
    t = pd.read_table(file_name, sep = ',')
    t = t.rename(columns=lambda x: x.strip())
    print t
    for gene, contig, plasmid, start, position in izip(t["qseqid"], t["contig"], t["sseqid"], t["sstart"], t["position"]):
        
        loc_map[str(gene)]={'gene':gene, 'contig':contig, 'plasmid':plasmid, 'bp_start': start, 'position': position}
    
    # Returning a tuple of the data frame and a dictionary where keys are the genes
    return (t, loc_map)

# Test
df_map, dict_map = (gene_locations())
print dict_map['group10']

def all_versus_corrected_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/',
                                plasmid_type = ['unitig_5']):

    nod_genes = OrderedDict([(4144, 'nodX'), (4143, 'nodN'), (4142, 'nodM'), (4141, 'nodL'), (4140, 'nodE'), (4139, 'nodF'), (4138, 'nodD'), 
    (4137, 'nodA'), (4136, 'nodC'), (4135, 'nodI'), (4134, 'nodJ'), (4129, 'nifB'), (4128, 'nifA'), (4127, 'fixX'), (4126, 'fixC'), (4125, 'fixB'), 
    (4124, 'fixA'), (4123, 'nifH'), (4122, 'nifD'), (4121, 'nifK'), (4120, 'nifE'), (2448, 'rpoB'), (2140, 'recA')])
    os.chdir(directory)

    #This is the gene SNPs matrix
    genes = []
    for f in glob.glob('*.npz'):
        with np.load(directory + f) as data:
            genes.append((f, data["matrix"], data["strains"]))

    r_scores = []
    gene_name = []
    original_name = []
    n_snps = []
    n_members = []
    origin = [] 

    print 'Number of genes being analized is: %d' % len(genes)

    # Initializing the LD table 

    # Subseting by a given condition (plasmid type)
    df_locations, dict_locations = (gene_locations())
    #sub_df = df_locations.loc[(df_locations.qseqid == plasmid_type)]

    LD_matrix = np.empty((len(df_locations['position']), len(df_locations['position'])))
    LD_matrix[:] = np.NAN
    LD_matrix = pd.DataFrame(LD_matrix, index = df_locations['qseqid'], columns = df_locations['qseqid'])
    
    print len(df_locations['position'])
    for index1 in range(0, len(genes)):
        for index2 in range(index1, len(genes)):
            
            gene1 = genes[index1]
            gene2 = genes[index2]
            
            name1, snps1, strains_1 = (gene1)
            name2, snps2, strains_2 = (gene2)

            # Find the plasmid origin of the gene:
            missing_location = 0 
            try:
                origin1 = dict_locations[name1[:-4]]['plasmid']
                origin2 = dict_locations[name2[:-4]]['plasmid']
            except KeyError:
                missing_location += 1
                continue

            if origin1 in plasmid_type and origin2 in plasmid_type:
            
                # Construct GRM for a singular gene
                grm_1 = np.divide(np.dot(snps1, snps1.T), snps_1.shape[1])

                flat_grm_1 = grm_1.flatten()
                norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
                norm_flat_grm1 = flat_grm_1
                norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))
            
                # Gene 2
                # Construct GRM for a singular gene
                grm_2 = np.divide(np.dot(snps2, snps2.T), snps2.shape[1])

                flat_grm_2 = grm_2.flatten()
                norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
                norm_flat_grm2 = flat_grm_2
                norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))
        
                corr = pearsonr(norm_flat_grm1, norm_flat_grm2)[0]
                
                print str(name1[:-4])
                print str(name2[:-4]) 
                print corr
                LD_matrix.set_value(name1[:-4], name2[:-4], corr)

                #print corr
                #print LD_matrix

    print ' Number of genes that have not a location assigned for the plasmid %s is: %d' % (plasmid_type, missing_location)

    LD_matrix.to_csv('LD_table_test.csv', header = True)
    return LD_matrix

print all_versus_corrected_genes()