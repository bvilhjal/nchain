"""
Code for analysing LD of rhizobium
"""

import scipy as sp
import h5py 
import os

nt_map = {'A':1, 'C':2, 'G':3, 'T':4, '-':5}

def gen_genotype_hdf5file(out_hdf5_file ='/faststorage/project/NChain/rhizobium/ld/snps2.hdf5', 
                          snps_directory='/faststorage/project/NChain/rhizobium/ld/snps/',
                          fna_files_directory='/faststorage/project/NChain/rhizobium/ld/group_alns/',
                          snps_wo_struct_directory='/faststorage/project/NChain/rhizobium/ld/snps_no_structure/'):

    h5f = h5py.File(out_hdf5_file)
    snps_g = h5f.create_group('snps')
    snps_wostruct_g = h5f.create_group('snps_wo_struct')
    align_g = h5f.create_group('alignments')
    
    print "Parsing alignments"
    aln_files = os.listdir(fna_files_directory)
    for i, a_f in enumerate(aln_files):
        if i%100==0:
            print i
        l = a_f.split('.')
        group_num = int(l[0][5:])
        aln_dict = parse_fasta_file(fna_files_directory+'/'+a_f)
        g = align_g.create_group('%d'%group_num)
        g.create_dataset('strains', data=aln_dict['iids'])
        g.create_dataset('sequences', data=aln_dict['sequences'], compression='lzf')
        g.create_dataset('psequences', data=aln_dict['psequences'], compression='lzf')
        g.create_dataset('nsequences', data=sp.array(aln_dict['nsequences'], dtype='int8'))

    print "Raw SNPs files"
    snps_files = os.listdir(snps_directory)
    for i, snps_f in enumerate(snps_files):
        if i%100==0:
            print i
        l = snps_f.split('.')
        group_num = int(l[0][5:])
        snps_data = sp.load(snps_directory+'/'+snps_f)
        g = snps_g.create_group('%d'%group_num)
        g.create_dataset('snps', data=snps_data['matrix'], compression='lzf')
        g.create_dataset('alignment_length',data=snps_data['alignment_length'])
        g.create_dataset('minor_frequencies',data=snps_data['minor_frequencies'])
        g.create_dataset('strains',data=[a.encode('utf8') for a in snps_data['strains']])
        g.create_dataset('positions',data=snps_data['positions'])
        h5f.flush()

    print "Now SNPs wo structure"
    snps_files = os.listdir(snps_wo_struct_directory)
    for i, snps_f in enumerate(snps_files):
        if i%100==0:
            print i
        l = snps_f.split('.')
        group_num = int(l[0][5:])
        snps_data = sp.load(snps_wo_struct_directory+'/'+snps_f)
        g = snps_wostruct_g.create_group('%d'%group_num)
        g.create_dataset('snps', data=snps_data['matrix'], compression='lzf')        
        g.create_dataset('alignment_length',data=snps_data['alignment_length'])
        g.create_dataset('minor_frequencies',data=snps_data['minor_frequencies'])
        g.create_dataset('strains',data=[a.encode('utf8') for a in snps_data['strains']])
        g.create_dataset('positions',data=snps_data['positions'])

    
    h5f.close()
    
    
    
def parse_blosum62(blosum62_file):
    from itertools import izip
    with open(blosum62_file) as f:
        blosum62_matrix = sp.zeros((24,24)) 
        blosum62_dict = {}
        mat_i = 0 
        for line in f:
            if line[0]=='#':
                continue
            elif line[0]==' ':
                aas = line.split()
                for aa in aas:
                    blosum62_dict[aa]={}
            else:
                l = line.split()
                aa1 = l[0]
                scores = map(float,l[1:])
                blosum62_matrix[mat_i]=sp.array(scores)
                for aa2, score in izip(aas,scores):
                    blosum62_dict[aa1][aa2]=score
                mat_i += 1
                
    return blosum62_matrix, blosum62_dict
                
            


def translate_dna(sequence):

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    '---':'-'
    }
    proteinsequence = ''
#     start = sequence.find('ATG')
#     sequencestart = sequence[int(start):]
#     stop = sequencestart.find('TAA')
#     cds = str(sequencestart[:int(stop)+3])
    cds=sequence
    
    for i in range(0,len(cds),3):
        proteinsequence += codontable[cds[i:i+3]]
    return proteinsequence

    
    
def parse_fasta_file(filename):
    
    header = None
    sequence = ''
    data_dict = {'iids':[], 'sequences':[], 'psequences':[], 'nsequences':[]}
    with open(filename) as f:
        for line in f:
            if line[0] == ">":
                if header is not None:
                    assert len(sequence)%3==0, 'Sequence length is not divisible by 3.'
                    l = header.split('|')
                    data_dict['iids'].append(l[2])
                    data_dict['sequences'].append(sequence)
                    psequence = translate_dna(sequence)
                    data_dict['psequences'].append(psequence)
                    nsequence = map(lambda x: nt_map[x],sequence)
                    data_dict['nsequences'].append(nsequence)
                header = line.strip()
                sequence = ''
            else:
                sequence += line.strip()
        
    return data_dict


    
def check_variants(gt_hdf5_file='snps2.hdf5'):
    # Iterate across all SNPs.
    # For each SNP 
    #    Check if it results in an AA change.
    #    If so quantify it.
    #    Count SNPs with more than 2 variants.
    #     
    pass
    
def call_variants(gt_hdf5_file='snps2.hdf5'):
    """
    Generate a new set of SNPs to look at.
    
    For all nts:
        if it is a SNP
            count # of variants. 
            check AA changes
            quantify AA change severity    
    """
    pass
    
    
    