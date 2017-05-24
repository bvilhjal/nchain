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
from Bio.Blast.NCBIXML import parse
import h5py
import glob
import pandas as pd
import numpy as np
import re

if argv[1] == 'mac':
    snps_file='/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5'
    out_dir='/Users/PM/Desktop/PHD_incomplete/Methods/Intergenic_LD/corrected_snps_test/'
    figure_dir='/Users/PM/Desktop/PHD_incomplete/nchain/Figures/' 
    alignments_dir = '/Users/PM/Desktop/PHD_incomplete/nchain/'
    database_file = '/Users/PM/Desktop/PHD_incomplete/Methods/Dot_plots/Polished_assemblies/'

if argv[1] == 'cluster':
    snps_file='/faststorage/project/NChain/rhizobium/ld/new_snps.hdf5'
    out_dir = '/faststorage/project/NChain/rhizobium/intergenic_LD/corrected_snps_test'
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

	# Blasting the gene sequence against its belonging scaffold

    query_name = 'query_' + genome + '.fasta'
    db_name = database_dir + genome + '.fasta'
    blastx_cline = NcbiblastnCommandline(query=str(query_name), db=str(db_name), evalue=0.1, outfmt = 6)
    out, err = blastx_cline()
    list_out = re.split('\n|\t', out)
    print out
    print list_out
    print len(list_out)
    del list_out[-1]
    blast_df = pd.DataFrame(np.array(list_out).reshape(len(list_out)/12,12), columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore'])

    # Changing the type of the data frame:
    blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']] = blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']].apply(pd.to_numeric)

    # Saving the blast results in csv file:
    name = 'blast_results' + genome
    blast_df.to_csv(name, sep = '\t')
    
    print blast_df.dtypes
    # Sorting the Dataframe by first the chromids and second by start (ascending True for both cases)
    blast_df_sorted = blast_df.sort_values(by = ['sseqid','start'], ascending= [True, True])
    print blast_df_sorted
blasting()




