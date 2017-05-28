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
from simple_intergenic_ld import correlation_plot
from sys import argv


# Setting the directory of the data depending on what computer I am  

if argv[1] == 'windows':
    file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/gene_locations_sorted.csv'
    
if argv[1] == 'mac':
    file_name ='/Users/PM/Desktop/PHD_incomplete/nchain/gene_locations_sorted.csv' 

if argv[1] == 'cluster':
    file_name = '/faststorage/project/NChain/rhizobium/intergenic_LD/'

def gene_locations(file_name = file_name):
    '''Adding the location of the genes based on a blast of the pacbio genome SM158'''
    
    loc_map = {}
    t = pd.read_table(file_name, sep = ',')
    t = t.rename(columns=lambda x: x.strip())
    for gene, contig, plasmid, start, position in izip(t["qseqid"], t["contig"], t["sseqid"], t["sstart"], t["position"]):
        
        loc_map[str(gene)]={'gene':gene, 'contig':contig, 'plasmid':plasmid, 'bp_start': start, 'position': position}
    
    # Returning a tuple of the data frame and a dictionary where keys are the genes
    return (t, loc_map)

# Test
df_map, dict_map = (gene_locations())
print dict_map['group10']

def all_versus_corrected_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps_test/',
                                plasmid_type = ['unitig_5'], formats = '.png'):

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

    # Subsetting by a given condition (plasmid type)
    df_locations, dict_locations = (gene_locations())

    sub_df = df_locations.loc[(df_locations['sseqid'] == plasmid_type[0])]

    LD_matrix = np.zeros((len(sub_df['position']), len(sub_df['position'])))
    LD_matrix = pd.DataFrame(LD_matrix, index = sub_df['qseqid'], columns = sub_df['qseqid'])
    
    print 'Number of genes in the data base %d' % len(df_locations['position'])
    missing_location = 0 
    for gene1 in genes:
        for gene2 in genes:
            
            #gene1 = genes[index1]
            #gene2 = genes[index2]
            
            name1, snps1, strains_1 = (gene1)
            name2, snps2, strains_2 = (gene2)

            # Find the plasmid origin of the gene:
            try:
                origin1 = dict_locations[name1[:-4]]['plasmid']
                origin2 = dict_locations[name2[:-4]]['plasmid']
            except KeyError:
                missing_location += 1
                continue

            if origin1 in plasmid_type and origin2 in plasmid_type:
            
                # This works only if we assume that strains_2 and strains_1 are ordered beforehand.  Are they? They are.
                strains_mask_1 = np.in1d(strains_1, strains_2, assume_unique=True)
                fitered_strains_1 = strains_1[strains_mask_1]
                strains_mask_2 = np.in1d(strains_2, fitered_strains_1, assume_unique=True)
                #fitered_strains_2 = strains_2[strain_mask_2]

                # Construct GRM for a singular gene
                total_snps_1 = snps1[strains_mask_1, :]
                grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])

                flat_grm_1 = grm_1.flatten()
                norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
                norm_flat_grm1 = flat_grm_1
                norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))
            
                # Gene 2
                # Construct GRM for a singular gene
                total_snps_2 = snps2[strains_mask_2, :]
                grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])

                flat_grm_2 = grm_2.flatten()
                norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
                norm_flat_grm2 = flat_grm_2
                norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))
        
                corr = pearsonr(norm_flat_grm1, norm_flat_grm2)[0]
                
                LD_matrix[str(name1[:-4])][str(name2[:-4])] = corr

                #print corr
                #print LD_matrix

    print 'Number of genes that have not a location assigned for the plasmid %s is: %d' % (plasmid_type, missing_location)

    # Removing rows with empty values
    LD_matrix = LD_matrix[(LD_matrix.T != 0).any()]

    # Removing columns with empty values
    LD_matrix = LD_matrix.loc[:, (LD_matrix != 0).any(axis = 0)]

    LD_matrix.to_csv(plasmid_type[0] + '.csv', header = True)

    correlation_plot(LD_matrix, wrt=False, fig_name = plasmid_type[0] + formats, show = False)
    return LD_matrix

print all_versus_corrected_genes(plasmid_type = ['unitig_0'])
#print all_versus_corrected_genes(plasmid_type = ['unitig_3'])
#print all_versus_corrected_genes(plasmid_type = ['unitig_2'])
#print all_versus_corrected_genes(plasmid_type = ['unitig_1'])
#print all_versus_corrected_genes(plasmid_type = ['unitig_0'])
