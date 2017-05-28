# Pipeline for the LD between different genomes

# 0 Create the query data file with the 5000 genes
# 1 Choose the Pacbio data
# 2 Find the position of the genes in the specific chosen genome
# 3 Calculate structure LD for each plasmid and between plasmids
# 4 Find LD blocks with very high correlation (in average aove 0.6)
# 5 Find functionality of those blocks


# 1 Chose Pacbio Data
# In the folder polished assemblies we have the needed files of each PacBio data to make the blast against the 5000 genes
from sys import argv 
from Bio.Blast.Applications import NcbiblastnCommandline
import h5py
import glob
import pandas as pd
import numpy as np
import scipy as sp
import re
from simple_intergenic_ld import correlation_plot
from sys import argv 
from itertools import izip
from collections import OrderedDict
from scipy.stats.stats import pearsonr
import os

# Setting the directory of the data depending on what computer I am  

if argv[1] == 'windows':
    snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'
    out_dir='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps_test/'
    figure_dir='C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Figures/'
    
if argv[1] == 'mac':
    snps_file='/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5'
    out_dir='/Users/PM/Desktop/PHD_incomplete/Methods/Intergenic_LD/corrected_snps_test/'
    figure_dir='/Users/PM/Desktop/PHD_incomplete/nchain/Figures/' 
    alignments_dir = '/Users/PM/Desktop/PHD_incomplete/nchain/'
    database_file = '/Users/PM/Desktop/PHD_incomplete/Methods/Dot_plots/Polished_assemblies/'

if argv[1] == 'cluster':
    snps_file='/faststorage/project/NChain/rhizobium/ld/new_snps.hdf5'
    out_dir = '/faststorage/project/NChain/rhizobium/intergenic_LD/corrected_snps_test/'
    figure_dir='/faststorage/project/NChain/rhizobium/intergenic_LD/figures' 
    alignments_dir = '/faststorage/project/NChain/rhizobium/Gene_group_Ontology/Blast_group_genes/'
    database_file = '/faststorage/project/NChain/rhizobium/pacbio_reads/'


pacbio_dict = {'SM170C_S17':['3405-170C', 'gsC'],
            'SM158_S16': ['3385-158', 'gsC'],
            'SM152B_S15': ['3370-152B', 'gsA'],
            'SM149A_S14': ['3362-149A', 'gsE'],
            'SM147A_S13': ['3356-147A','gsC'],
            'SM41_S12': ['3239-41','gsC'],
            'SM3_S10': ['3206-3', 'gsB'] }   

genome = argv[2]
# Creating the query data file with the 5000 genes

def parse_fasta(filename, genome):
    file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
    file_separe = file.split('>') #spliting each entry by the > 
    file_separe.remove('')
    parse_dict = {}
    header = []
    for entry in file_separe:
        seq = entry.splitlines()
        header = seq[0] #these are the first elements of the list 
        
        genome_id = seq[0].split('|')[1]

        if pacbio_dict[genome][0] == genome_id:
            seq = ''.join(seq[1:]) #joining the sequences 
            parse_dict[header] = seq
    return parse_dict

def create_query(snp_file, align_file, genome = genome):

    '''Download the alignments of each gene group and make a big query fasta file'''

    gene_parse = list()
    fasta_name = 'query_'+ genome 

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    for gene in gene_groups:
        file_name = align_file + 'group' + gene + '.fna'
        gene_parse.append(parse_fasta(filename = file_name, genome = genome))
    print gene_parse
    with open('%s.fasta' % fasta_name, 'w') as f:
        for gene in gene_parse:
            for names, seq in gene.items():
                f.write('>{}\n'.format(names))
                f.write('{}\n\n'.format(seq))
    f.close()
    #return gene_parse

#print create_query(snp_file = snps_file, align_file = alignments_dir, genome = genome)

def blasting(genome = genome, evalue = 0.001, database_dir = database_file):
    
    #create_query(snp_file = snps_file, align_file = alignments_dir, genome = genome)
    # Blasting the gene sequence against its belonging scaffold

    query_name = 'query_' + genome + '.fasta'
    db_name = database_dir + genome + '.fasta'
    blastx_cline = NcbiblastnCommandline(query=str(query_name), db=str(db_name), evalue=0.1, outfmt = 6, max_target_seqs = 1)
    out, err = blastx_cline()
    list_out = re.split('\n|\t', out)
    print out
    print list_out
    print len(list_out)
    del list_out[-1]
    blast_df = pd.DataFrame(np.array(list_out).reshape(len(list_out)/12,12), columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore'])

    # Changing the type of the data frame:
    blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']] = blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']].apply(pd.to_numeric)

    # Deleting duplicates by a given column (evalue)
    e_maxes = blast_df.groupby(['qseqid']).evalue.transform(max)
    blast_df = blast_df[blast_df.evalue == e_maxes]

    # Saving the blast results in csv file:
    name = 'blast_results_test' + genome + '.csv'
    blast_df.to_csv(name, sep = '\t')
    
    print blast_df.dtypes
    # Sorting the Dataframe by first the chromids and second by start (ascending True for both cases)
    blast_df_sorted = blast_df.sort_values(by = ['sseqid','start'], ascending= [True, True])

    # Spliting the first column:
    blast_df_sorted['gene'], blast_df_sorted['strain'], blast_df_sorted['contig'] = zip(*blast_df_sorted['qseqid'].map(lambda x: x.split('|')))
    
    print blast_df_sorted
    file_name = 'sorted_genes_testing' + genome +'.csv'
    blast_df_sorted.to_csv(file_name, sep=',')
    return file_name

def gene_locations(file_name):
    '''Adding the location of the genes based on a blast of the pacbio genome SM158'''

    loc_map = {}
    t = pd.read_table(file_name, sep = ',')
   # print t
    t = t.rename(columns=lambda x: x.strip())
    #print t
    for gene, contig, plasmid, start in izip(t["gene"], t["contig"], t["sseqid"], t["start"]):
        loc_map[str(gene)]={'gene':gene, 'contig':contig, 'plasmid':plasmid, 'bp_start': start}
    
    # Returning a tuple of the data frame and a dictionary where keys are the genes
    return (t, loc_map)

def all_versus_corrected_genes(directory = out_dir,
                                plasmid_type = [argv[3]], formats = '.png', genome = genome):

    nod_genes = OrderedDict([(4144, 'nodX'), (4143, 'nodN'), (4142, 'nodM'), (4141, 'nodL'), (4140, 'nodE'), (4139, 'nodF'), (4138, 'nodD'), 
    (4137, 'nodA'), (4136, 'nodC'), (4135, 'nodI'), (4134, 'nodJ'), (4129, 'nifB'), (4128, 'nifA'), (4127, 'fixX'), (4126, 'fixC'), (4125, 'fixB'), 
    (4124, 'fixA'), (4123, 'nifH'), (4122, 'nifD'), (4121, 'nifK'), (4120, 'nifE'), (2448, 'rpoB'), (2140, 'recA')])
    os.chdir(directory)
    
    file_saved = blasting()
    df_locations, dict_locations = (gene_locations(file_saved))
    
    gene_keys = dict_locations.keys()
    #This is the gene SNPs matrix
    genes = []
    for f in gene_keys:
 
        name = directory + f + '.npz'
        with np.load(name) as data:
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
    sub_df = df_locations.loc[(df_locations['sseqid'] == plasmid_type[0])]

    print sub_df

    LD_matrix = np.zeros((len(sub_df['contig']), len(sub_df['contig'])))
    LD_matrix = pd.DataFrame(LD_matrix, index = sub_df['gene'], columns = sub_df['gene'])
    index_names = list(LD_matrix.index.values)

    print 'Number of genes in the data base %d' % len(df_locations['gene'])
    missing_location = 0 
    for gene1 in genes:
        for gene2 in genes:
            
            name1, snps1, strains_1 = (gene1)
            name2, snps2, strains_2 = (gene2)
            
            # Find the plasmid origin of the gene:
            try:
                origin1 = dict_locations[name1]['plasmid']
                origin2 = dict_locations[name2]['plasmid']
            except KeyError:
                missing_location += 1
                continue
    
            if origin1 in plasmid_type and origin2 in plasmid_type:
                
                print origin1 
                print origin2
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
                
                l = index_names.index(str(name1))
                c = index_names.index(str(name2))
                
                LD_matrix.ix[l,c] = corr
    print 'Number of genes that have not a location assigned for the plasmid %s is: %d' % (plasmid_type, missing_location)
    
    # Removing rows with empty values
    LD_matrix = LD_matrix[(LD_matrix.T != 0).any()]
    
    # Removing columns with empty values
    LD_matrix = LD_matrix.loc[:, (LD_matrix != 0).any(axis = 0)]

    LD_matrix.to_csv(plasmid_type[0] + '.csv', header = True)

    correlation_plot(LD_matrix, wrt=False, fig_name = plasmid_type[0] + formats, show = False)
    return LD_matrix

print all_versus_corrected_genes(plasmid_type = [arvg[3]])
