"""
Code for analysing LD of rhizobium
"""

import scipy as sp
import h5py 
import os
import matplotlib
matplotlib.use('Agg')
import pylab

nt_map = {'A':1, 'C':2, 'G':3, 'T':4, '-':5, 'N':6}
nt_decode_map = {1:'A', 2:'C', 3:'G', 4:'T', 5:'-', 6:'N'}

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
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    '---':'X'
    }



def get_codon_syn_map():
    all_codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 
                  'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 
                  'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 
                  'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 
                  'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 
                  'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 
                  'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 
                  'TTT']
    
    all_bases = {'A','T','C','G'}

    ret_dict = {}
    
    for codon in all_codons:
        syn_list = []
        for i in range(3):
            base = codon[i]
            other_bases = all_bases - {base}
            syn = 0
            for new_base in other_bases:
                new_codon = codon[:i] + new_base + codon[i + 1:]
                syn += int(codontable[codon]==codontable[new_codon])
            syn_list.append(syn/3.0)
        ret_dict[codon]=sp.sum(syn_list)
    return ret_dict


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

    proteinsequence = ''
#     start = sequence.find('ATG')
#     sequencestart = sequence[int(start):]
#     stop = sequencestart.find('TAA')
#     cds = str(sequencestart[:int(stop)+3])
    cds=sequence
    
    for i in range(0,len(cds),3):
        codon = cds[i:i+3]

        proteinsequence += codontable.get(codon,'X')
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
    
    
    
def call_variants(gt_hdf5_file='snps2.hdf5', out_file='new_snps.hdf5', min_num_strains=100, 
                  blosum62_file='/home/bjarni/NChain/faststorage/rhizobium/ld/blosum62.txt'):
    """
    Generate a new set of SNPs to look at.
    
    For all nts:
        if it is a SNP
            count # of variants. 
            check AA changes
            quantify AA change severity    
    
    """
    from itertools import izip
    blosum62_matrix, blosum62_dict = parse_blosum62(blosum62_file)
    codon_syn_map = get_codon_syn_map()
    h5f = h5py.File(gt_hdf5_file)
    ag = h5f['alignments']
    oh5f = h5py.File(out_file)
    gene_groups = sorted(ag.keys())
    num_parsed_genes = 0
    for gg in gene_groups:
        g = ag[gg]
        if len(g['strains'])>min_num_strains:
            
            #1. Filter indel/bad rows
            nt_mat = g['nsequences'][...]
            num_vars = sp.apply_along_axis(lambda x: len(sp.unique(x)), 0, nt_mat)
            nt_mat = sp.transpose(nt_mat)
            bad_rows_filter = num_vars<5
            if sp.sum(bad_rows_filter)>0:
                raw_snps = nt_mat[bad_rows_filter]
                
                #2. Filter non-variable rows
                ok_num_vars = num_vars[bad_rows_filter]
                var_filter = ok_num_vars>1                
                num_raw_snps = sp.sum(var_filter)
                if num_raw_snps>0:
                    print 'Working on gene group: %s'%gg
                    
                    M,N = nt_mat.shape
                    var_positions = sp.arange(M)[bad_rows_filter]
                    all_snps = raw_snps[var_filter]
                    all_snp_positions = var_positions[var_filter]
                    
                    #3. Identify good SNPs (dimorphic SNPs)
                    good_snp_filter = ok_num_vars==2
                    ok_snps = raw_snps[good_snp_filter]
                    snp_positions = var_positions[good_snp_filter]
                    assert len(ok_snps)==len(snp_positions), 'A bug detected!'
                    
                    #4. Call good SNPs
                    sequence = g['sequences'][0]
                    snps = []
                    nts = []
                    
                    codon_snps = []
                    codon_snp_positions = []
                    codons = []
                    aacids = []
                    is_synonimous_snp =  []
                    blosum62_scores = []
                    tot_num_syn_sites = 0
                    tot_num_non_syn_sites = 0
                    for ok_snp, snp_pos in izip(ok_snps, snp_positions):                    
                        mean_snp = sp.mean(ok_snp)
                        snp = sp.zeros(N)
                        snp[ok_snp>mean_snp]=1
                        snps.append(snp)
                        
                        #Get nucleotides 
                        nt0 = nt_decode_map[ok_snp.min()]
                        nt1 = nt_decode_map[ok_snp.max()]
                        nts.append([nt0,nt1])
                        
                        #5. Check codon position
                        codon_pos = snp_pos%3
                        
                        #6. Identify non-ambiguous codon changes.
                        #Check if there is a preceding/succeeding SNP within the codon.
                        if codon_pos==0:
                            if not (bad_rows_filter[snp_pos+1] and bad_rows_filter[snp_pos+2]):
                                continue
                            if not(num_vars[snp_pos+1]==1 and num_vars[snp_pos+2]==1):
                                continue
                            cdseq12 = sequence[snp_pos+1:snp_pos+3]
                            codon0 = nt0+cdseq12
                            codon1 = nt1+cdseq12
                            
                        elif codon_pos==1:
                            if not (bad_rows_filter[snp_pos-1] and bad_rows_filter[snp_pos+1]):
                                continue
                            if not(num_vars[snp_pos-1]==1 and num_vars[snp_pos+1]==1):
                                continue
                            cdseq0 = sequence[snp_pos-1]
                            cdseq2 = sequence[snp_pos+1]
                            codon0 = cdseq0+nt0+cdseq2
                            codon1 = cdseq0+nt1+cdseq2
                        
                        elif codon_pos==2:
                            if not (bad_rows_filter[snp_pos-1] and bad_rows_filter[snp_pos-2]):
                                continue
                            if not(num_vars[snp_pos-1]==1 and num_vars[snp_pos-2]==1):
                                continue
                            cdseq01 = sequence[snp_pos-2:snp_pos]
                            codon0 = cdseq01+nt0
                            codon1 = cdseq01+nt1
                        
                        assert codon0!=codon1, 'Codons are identical?'
                        
                        #This appears to be a unique codon change with a dimorphic SNP.
                        codons.append([codon0,codon1])
                        freq = sp.mean(snp,0)
                        
                        #7. Check non-synonimous/synonimous
                        num_syn_sites = freq*codon_syn_map[codon0]+(1-freq)*codon_syn_map[codon1]
                        num_non_syn_sites = 3-num_syn_sites
                        tot_num_syn_sites += num_syn_sites
                        tot_num_non_syn_sites += num_non_syn_sites

                        aa0 = codontable[codon0]
                        aa1 = codontable[codon1]
                        aacids.append([aa0,aa1])
                        is_synon = aa0==aa1
                        is_synonimous_snp.append(is_synon)
                        
                        
                        #Get BLOUSUM62 score
                        blosum62_score = blosum62_dict[aa0][aa1]
                        blosum62_scores.append(blosum62_score)

                        codon_snps.append(snp)
                        codon_snp_positions.append(snp_pos)
                    
                    #Normalize SNPs
                    freqs = sp.mean(snps,0)
                    norm_snps = (snps-freqs)/sp.sqrt(freqs*(1-freqs))
                    
                    #Calculate dn/ds ratios
                    num_syn_subt = sp.sum(is_synonimous_snp)
                    num_non_syn_subt = len(is_synonimous_snp)-num_syn_subt
                    if num_non_syn_subt>0:
                        dn_ds_ratio = (num_non_syn_subt/tot_num_non_syn_sites)/(num_syn_subt/tot_num_syn_sites)
                    else:
                        dn_ds_ratio=-1

                    
                    #Store everything to a HDF5 file
                    og = oh5f.create_group(gg)   
                    og.create_dataset('var_positions', data=var_positions)
                    og.create_dataset('num_vars', data=num_vars)
                    og.create_dataset('raw_snps', data=sp.array(all_snps,dtype='int8'), compression='lzf')
                    og.create_dataset('raw_snp_positions', data=all_snp_positions)
                    og.create_dataset('snps', data=sp.array(snps,dtype='int8'), compression='lzf')
                    og.create_dataset('norm_snps', data=sp.array(norm_snps,dtype='single'), compression='lzf')
                    og.create_dataset('snp_positions', data=snp_positions)
                    og.create_dataset('codon_snps', data=codon_snps)
                    og.create_dataset('codon_snp_positions', data=codon_snp_positions)
                    og.create_dataset('blosum62_scores', data=blosum62_scores)
                    og.create_dataset('aacids', data=sp.array(aacids))
                    og.create_dataset('nts', data=sp.array(nts))
                    og.create_dataset('codons', data=sp.array(codons))
                    og.create_dataset('num_syn_sites', data=tot_num_syn_sites)
                    og.create_dataset('num_non_syn_sites', data=tot_num_non_syn_sites)
                    og.create_dataset('dn_ds_ratio', data=dn_ds_ratio)
                    oh5f.flush()
                    num_parsed_genes +=1

    print 'Parsed %d'%num_parsed_genes
    
    
def summarize_nonsynonimous_snps(snps_hdf5_file, fig_dir):
    h5f = h5py.File(snps_hdf5_file)
    gene_groups = sorted(h5f.keys())
    dn_ds_ratios = []
    mean_blosum_62_scores = []
    for gg in gene_groups:
        g = h5f[gg]
        codon_snp_positions = g['codon_snp_positions'][...]
        if len(codon_snp_positions)>100:
            dn_ds_ratio = g['dn_ds_ratio'][...]
            if dn_ds_ratio==-1:
                dn_ds_ratios.append(0)
            else:
                dn_ds_ratios.append(1/dn_ds_ratio)

            blosum62_scores = sp.mean(g['blosum62_scores'][...])
            dn_ds_ratios.append(dn_ds_ratio)
            mean_blosum_62_scores.append(sp.mean(blosum62_scores))
    
    mean_blosum_62_scores = sp.nan_to_num(mean_blosum_62_scores)
    print 'Average dn/ds ration: %0.4f'%sp.mean(dn_ds_ratios)
    pylab.hist(dn_ds_ratios, bins=100)
    pylab.savefig(fig_dir+'/dn_ds_ratio.png')
        
    pylab.clf()
    pylab.hist(mean_blosum_62_scores)
    pylab.savefig(fig_dir+'/mean_blosum_62_scores.png')
       
    
    