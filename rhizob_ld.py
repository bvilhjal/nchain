"""
Code for analysing LD of rhizobium
"""

import scipy as sp
import h5py 
import os

nt_map = {'A':1, 'C':2, 'G':3, 'T':4, '-':5}

def gen_genotype_hdf5file(out_hdf5_file ='/faststorage/project/NChain/rhizobium/ld/snps.hdf5', 
                          snps_directory='/faststorage/project/NChain/rhizobium/ld/snps/',
                          fna_files_directory='/faststorage/project/NChain/rhizobium/ld/group_alns/',
                          snps_wo_struct_directory='/faststorage/project/NChain/rhizobium/ld/snps_no_structure/'):

    h5f = h5py.File(out_hdf5_file)
    g = h5f.create_group('gene_groups')
    snps_g = g.create_group('snps')
    snps_wostruct_g = g.create_group('snps_wo_struct')
    align_g = g.create_group('alignments')
    
    print "Raw SNPs files"
    snps_files = os.listdir(snps_directory)
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
    }
    proteinsequence = ''
    start = sequence.find('ATG')
    sequencestart = sequence[int(start):]
    stop = sequencestart.find('TAA')
    cds = str(sequencestart[:int(stop)+3])

    for n in range(0,len(cds),3):
        if cds[n:n+3] in codontable == True:
            proteinsequence += codontable[cds[n:n+3]]
            print proteinsequence
        sequence = ''


    
    
def parse_fasta_file(filename):
    
    header = None
    sequence = None
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
                    data_dict['nsequence'].append(nsequence)
                header = line.strip()
                sequence = ''
            else:
                sequence += line.strip()
        
    return data_dict
    
    
    
    