"""
Code for analysing LD of rhizobium
"""

import scipy as sp
import h5py 
import os
import matplotlib
matplotlib.use('Agg')
import pylab
import collections
import pandas as pd

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
        
        #0. Check if there is evidence for CNVs/paralogs?
        seq_ids = g['strains']
        strains_list = map(lambda x: x.split('-')[0], seq_ids)
        strains, strain_counts = sp.unique(strains_list, return_counts=True)
        if len(strains)<len(strains_list):
            print 'Evidence for paralogs/CNVs'
        elif len(seq_ids)>min_num_strains:
            strains = map(lambda x: x.split('-')[0], seq_ids)
                        
            #1. Filter indel/bad rows
            nt_mat = g['nsequences'][...]
            num_vars = sp.apply_along_axis(lambda x: len(sp.unique(x)), 0, nt_mat)
            nt_mat = sp.transpose(nt_mat)
            bad_rows_filter = num_vars<5
            if sp.sum(bad_rows_filter)>0:
                raw_snps = nt_mat[bad_rows_filter]
                
                #Calculate nucleotide diversity
                M,N = raw_snps.shape
                diversity = 0.0
                for i in range(N-1):
                    for j in range(i+1,N):
                        diversity+=sp.sum(raw_snps[:,i]!=raw_snps[:,j])
                
                diversity = diversity/len(raw_snps)
                diversity = 2*diversity/(N*(N-1))
                
                #2. Filter non-variable rows
                ok_num_vars = num_vars[bad_rows_filter]
                var_filter = ok_num_vars>1                
                num_raw_snps = sp.sum(var_filter)
                if num_raw_snps>0:
                    print 'Working on gene group: %s'%gg
                    
                    M,N = nt_mat.shape
                    non_gap_positions = sp.arange(M)[bad_rows_filter]
                    all_snps = raw_snps[var_filter]
                    all_snp_positions = non_gap_positions[var_filter]
                    
                    #3. Identify good SNPs (dimorphic SNPs)
                    good_snp_filter = ok_num_vars==2
                    ok_snps = raw_snps[good_snp_filter]
                    snp_positions = non_gap_positions[good_snp_filter]
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
                    norm_snps = sp.transpose(snps)
                    freqs = sp.mean(norm_snps,0)
                    norm_snps = (norm_snps-freqs)/sp.sqrt(freqs*(1-freqs))
                    norm_snps = sp.transpose(norm_snps)
                    
                    
                    norm_codon_snps = sp.transpose(codon_snps)
                    codon_snp_freqs = sp.mean(norm_codon_snps,0)
                    norm_codon_snps = (norm_codon_snps-codon_snp_freqs)/sp.sqrt(codon_snp_freqs*(1-codon_snp_freqs))
                    norm_codon_snps = sp.transpose(norm_codon_snps)

                    
                    #Calculate dn/ds ratios
                    num_syn_subt = sp.sum(is_synonimous_snp)
                    num_non_syn_subt = len(is_synonimous_snp)-num_syn_subt
                    if num_syn_subt>0:
                        dn_ds_ratio = (num_non_syn_subt/tot_num_non_syn_sites)/(num_syn_subt/tot_num_syn_sites)
                    else:
                        dn_ds_ratio=-1

                    
                    #Store everything to a HDF5 file
                    og = oh5f.create_group(gg)   
                    og.create_dataset('num_vars', data=num_vars)
                    og.create_dataset('raw_snps', data=sp.array(all_snps,dtype='int8'), compression='lzf')
                    og.create_dataset('raw_snp_positions', data=all_snp_positions)
                    og.create_dataset('snps', data=sp.array(snps,dtype='int8'), compression='lzf')
                    og.create_dataset('norm_snps', data=sp.array(norm_snps,dtype='single'), compression='lzf')
                    og.create_dataset('freqs', data=sp.array(freqs,dtype='single'))
                    og.create_dataset('snp_positions', data=snp_positions)
                    og.create_dataset('codon_snps', data=sp.array(codon_snps,dtype='single'), compression='lzf')
                    og.create_dataset('norm_codon_snps', data=sp.array(norm_codon_snps,dtype='single'), compression='lzf')
                    og.create_dataset('codon_snp_freqs', data=sp.array(codon_snp_freqs,dtype='single'))
                    og.create_dataset('is_synonimous_snp', data=is_synonimous_snp)
                    og.create_dataset('strains', data=strains)
                    og.create_dataset('codon_snp_positions', data=codon_snp_positions)
                    og.create_dataset('blosum62_scores', data=blosum62_scores)
                    og.create_dataset('aacids', data=sp.array(aacids))
                    og.create_dataset('nts', data=sp.array(nts))
                    og.create_dataset('codons', data=sp.array(codons))
                    og.create_dataset('num_syn_sites', data=tot_num_syn_sites)
                    og.create_dataset('num_non_syn_sites', data=tot_num_non_syn_sites)
                    og.create_dataset('dn_ds_ratio', data=dn_ds_ratio)
                    og.create_dataset('diversity', data=diversity)
                    oh5f.flush()
                    num_parsed_genes +=1
        else:
            print 'Too few strains..'
    print 'Parsed %d'%num_parsed_genes
    
    
    
def summarize_nonsynonimous_snps(snps_hdf5_file = '/project/NChain/faststorage/rhizobium/ld/called_snps.hdf5', 
                                 seq_file = '/project/NChain/faststorage/rhizobium/ld/snps2.hdf5',
                                 fig_dir = '/project/NChain/faststorage/rhizobium/ld/figures',
                                 geno_species='gsA'):
    h5f = h5py.File(snps_hdf5_file)
    sh5f = h5py.File(seq_file)
    gene_groups = sorted(h5f.keys())
    dn_ds_ratios = []
    mean_blosum_62_scores = []
    num_seg_sites = []
    pi_diversity = []
    for gg in gene_groups:
        g = h5f[gg]
        sg = sh5f['snps'][gg]
        codon_snp_positions = g['codon_snp_positions'][...]
        if len(codon_snp_positions)>100:
            dn_ds_ratio = g['dn_ds_ratio'][...]
            dn_ds_ratios.append(dn_ds_ratio)

            blosum62_scores = sp.mean(g['blosum62_scores'][...])
            mean_blosum_62_scores.append(sp.mean(blosum62_scores))
        
            raw_snp_positions = g['raw_snp_positions'][...]
            num_seg_sites_per_base = len(raw_snp_positions)/float(sg['alignment_length'][...])
            num_seg_sites.append(num_seg_sites_per_base)
            
            diversity = g['diversity'][...]
            pi_diversity.append(diversity)
    
    mean_blosum_62_scores = sp.nan_to_num(mean_blosum_62_scores)
    print 'Average dn/ds ration: %0.4f'%sp.mean(dn_ds_ratios)
    pylab.clf()
    pylab.hist(dn_ds_ratios, bins=100)
    pylab.title(r'$\frac{d_n}{d_s}$ (values below 1 suggest purifying selection.)')
    pylab.savefig(fig_dir+'/dn_ds_ratio.png')
        
#     pylab.clf()
#     pylab.hist(mean_blosum_62_scores, bins=100)
#     pylab.title('Average BLOSUM62 scores (values above 1 suggest purifying selection)')    
#     pylab.savefig(fig_dir+'/mean_blosum_62_scores.png')
    
    pylab.clf()
    pylab.hist(num_seg_sites, bins=100)
    pylab.title(r'Number of segregating sites per nucleotide ($S$)')    
    pylab.savefig(fig_dir+'/segregating_sites.png')
    
    pylab.clf()
    pylab.hist(pi_diversity, bins=100)
    pylab.title(r'Nucleotide diversity ($\pi$)')    
    pylab.savefig(fig_dir+'/nucleotide_diversity.png')
    
    pop_map, ct_array = parse_pop_map()
    unique_gs = sp.unique(ct_array)
    avg_gene_genosp_ld_dict = gene_genospecies_corr()
    print avg_gene_genosp_ld_dict
    #Now plot things per genospecies...
    dn_ds_ratios = {'all':[], 'nonsyn':[], 'syn':[]}
    num_seg_sites = {'all':[], 'nonsyn':[], 'syn':[]}
    pi_diversity = {'all':[], 'nonsyn':[], 'syn':[]}
    mean_r2s = {'all':[], 'nonsyn':[], 'syn':[]}
    
    for gg in gene_groups:
        g = h5f[gg]
        for snp_type in ['all','nonsyn','syn']:
            ok_genes = set(avg_gene_genosp_ld_dict[snp_type].keys())
            if gg in ok_genes:
                mean_r2s[snp_type].append(avg_gene_genosp_ld_dict[snp_type][gg][geno_species])
                dn_ds_ratio = g['dn_ds_ratio'][...]
                dn_ds_ratios[snp_type].append(dn_ds_ratio)
              
                
                raw_snp_positions = g['raw_snp_positions'][...]
                num_seg_sites_per_base = len(raw_snp_positions)/float(sg['alignment_length'][...])
                num_seg_sites[snp_type].append(num_seg_sites_per_base)
                 
                diversity = g['diversity'][...]
                pi_diversity[snp_type].append(diversity)

    
    for snp_type in ['all','nonsyn','syn']:
        pylab.clf()
        pylab.plot(mean_r2s[snp_type], dn_ds_ratios[snp_type], 'k.', alpha=0.3)
        pylab.xlabel('Mean r2 between SNPs within a gene and %s'%geno_species)
        pylab.ylabel(r'$\frac{K_a}{K_s}$')
        pylab.savefig(fig_dir+'/Ka_Ks_vs_%s_corr_%s.png'%(geno_species,snp_type))
        
        pylab.clf()
        pylab.plot(mean_r2s[snp_type], num_seg_sites[snp_type], 'k.', alpha=0.3)
        pylab.xlabel('Mean r2 between SNPs within a gene and %s'%geno_species)
        pylab.ylabel(r'Number of segregating sites per nucleotide ($S$)')
        pylab.savefig(fig_dir+'/Num_seg_sites_vs_%s_corr_%s.png'%(geno_species,snp_type))
        
        pylab.clf()
        pylab.plot(mean_r2s[snp_type], num_seg_sites[snp_type], 'k.', alpha=0.3)
        pylab.xlabel('Mean r2 between SNPs within a gene and %s'%geno_species)
        pylab.ylabel(r'Nucleotide diversity ($\pi$)')
        pylab.savefig(fig_dir+'/nt_diversity_vs_%s_corr_%s.png'%(geno_species,snp_type))

    
    
def gene_genospecies_corr(snps_hdf5_file = '/project/NChain/faststorage/rhizobium/ld/called_snps.hdf5',
                          min_maf = 0.15, min_num_snps = 20):
    from itertools import izip
    h5f = h5py.File(snps_hdf5_file)
    gene_groups = sorted(h5f.keys())

    pop_map, ct_array = parse_pop_map()
    unique_gs = sp.unique(ct_array)  
    avg_gene_genosp_ld_dict = {'all': {}, 'nonsyn': {}, 'syn': {}}
        
    for i, gg in enumerate(gene_groups):
        if i%100==0:
            print '%d: Gene %s'%(i,gg)  
        g = h5f[gg]

        #Filtering SNPs with small MAFs
        freqs = g['codon_snp_freqs'][...]
        mafs = sp.minimum(freqs,1-freqs)
        maf_filter = mafs>min_maf
        if sp.sum(maf_filter)>min_num_snps:            
            is_synonimous_snp = g['is_synonimous_snp'][...]
            is_nonsynonimous_snp = sp.negative(is_synonimous_snp)
            syn_snp_filter = is_synonimous_snp*maf_filter
            nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter

            if sp.sum(syn_snp_filter)>sp.sum(nonsyn_snp_filter):
                all_norm_snps = g['norm_codon_snps'][...]
                norm_snps = all_norm_snps[maf_filter]
                M,N = norm_snps.shape
                strains = g['strains']
                gs_list = sp.array([pop_map.get(strain,'NA') for strain in strains])
                
                d = {}
                for gs in unique_gs:          
                    gs_snp = sp.array(sp.in1d(gs_list,[gs]),dtype='single')
                    gs_snp = (gs_snp - sp.mean(gs_snp))/sp.std(gs_snp)
                    r_list = sp.dot(gs_snp,norm_snps.T)/float(N)
                    r2_list = r_list**2
                    assert M==len(r2_list), 'A bug detected.'
                    d[gs] = {'mean_r2': sp.mean(r2_list), 'num_snps': M, 'r2s': r2_list}
                avg_gene_genosp_ld_dict['all'][gg] = d
                
                nonsyn_snp_filter = is_synonimous_snp[maf_filter]
                M = sp.sum(nonsyn_snp_filter)
                if M>10:
                    d = {}
                    for gs in unique_gs:          
                        gs_snp = sp.array(sp.in1d(gs_list,[gs]),dtype='single')
                        gs_snp = (gs_snp - sp.mean(gs_snp))/sp.std(gs_snp)
                        r2_list = avg_gene_genosp_ld_dict['all'][gg][gs]['r2s'][nonsyn_snp_filter]
                        assert M==len(r2_list), 'A bug detected.'
                        d[gs] = {'mean_r2': sp.mean(r2_list), 'num_snps': M, 'r2s': r2_list}
                    avg_gene_genosp_ld_dict['nonsyn'][gg] = d
                    
                syn_snp_filter = sp.negative(nonsyn_snp_filter)
                M = sp.sum(syn_snp_filter)
                if sp.sum(syn_snp_filter)>10:
                    d = {}
                    for gs in unique_gs:          
                        gs_snp = sp.array(sp.in1d(gs_list,[gs]),dtype='single')
                        gs_snp = (gs_snp - sp.mean(gs_snp))/sp.std(gs_snp)
                        r2_list = avg_gene_genosp_ld_dict['all'][gg][gs]['r2s'][syn_snp_filter]
                        assert M==len(r2_list), 'A bug detected.'
                        d[gs] = {'mean_r2': sp.mean(r2_list), 'num_snps': M, 'r2s': r2_list}
                    avg_gene_genosp_ld_dict['syn'][gg] = d
    return avg_gene_genosp_ld_dict
    


def gen_ld_plots(snps_hdf5_file = '/project/NChain/faststorage/rhizobium/ld/called_snps.hdf5', 
                 max_dist=3000, min_maf=0.1, bin_size=60,
                 fig_dir = '/project/NChain/faststorage/rhizobium/ld', filter_pop=None):
    
    pop_map, ct_array = parse_pop_map()

    from itertools import izip
    h5f = h5py.File(snps_hdf5_file)
    gene_groups = sorted(h5f.keys())
    ld_dist_dict = {'all':{}, 'nonsyn':{}, 'syn':{}}
    distances = range(1,max_dist)
    
    for dist in distances:
        ld_dist_dict['all'][dist]={'r2_sum':0.0, 'snp_count':0.0}
        ld_dist_dict['nonsyn'][dist]={'r2_sum':0.0, 'snp_count':0.0}
        ld_dist_dict['syn'][dist]={'r2_sum':0.0, 'snp_count':0.0}
    
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
                        if pop_map[s]==filter_pop:
                            indiv_filter[s_i]=True
                    except:
                        continue
                if sp.sum(indiv_filter)<20:
                    continue
    
            
            if filter_pop is not None:
                codon_snps = g['codon_snps'][...]
                codon_snps = codon_snps[:,indiv_filter]
                norm_codon_snps = sp.transpose(codon_snps)
                freqs = sp.mean(norm_codon_snps,0)
                norm_codon_snps = (norm_codon_snps-freqs)/sp.sqrt(freqs*(1-freqs))
                norm_codon_snps = sp.transpose(norm_codon_snps)
                mafs = sp.minimum(freqs,1-freqs)
                maf_filter = mafs>min_maf
                if sp.sum(maf_filter)>1:
                    is_synonimous_snp = g['is_synonimous_snp'][...]
                    is_nonsynonimous_snp = sp.negative(is_synonimous_snp)
                    syn_snp_filter = is_synonimous_snp*maf_filter
                    nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter
    
                    if sp.sum(syn_snp_filter)>sp.sum(nonsyn_snp_filter):
                        all_norm_snps = norm_codon_snps
                        all_positions = g['codon_snp_positions'][...]
                        norm_snps = all_norm_snps[maf_filter]
                        positions = all_positions[maf_filter]
                        M,N = norm_snps.shape
                        
                        ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
                        assert M==len(positions), 'A bug detected.'
                        for i in range(M-1):
                            for j in range(i+1,M):
                                dist = positions[j] - positions[i]
                                if dist<max_dist:
                                    ld_dist_dict['all'][dist]['r2_sum']+=ld_mat[i,j]**2
                                    ld_dist_dict['all'][dist]['snp_count']+=1.0
        
                        
                        if sp.sum(nonsyn_snp_filter)>1:
                            norm_snps = all_norm_snps[nonsyn_snp_filter]
                            positions = all_positions[nonsyn_snp_filter]
                            M,N = norm_snps.shape
                            
                            ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
                            assert M==len(positions), 'A bug detected.'
                            for i in range(M-1):
                                for j in range(i+1,M):
                                    dist = positions[j] - positions[i]
                                    if dist<max_dist:
                                        ld_dist_dict['nonsyn'][dist]['r2_sum']+=ld_mat[i,j]**2
                                        ld_dist_dict['nonsyn'][dist]['snp_count']+=1.0
                                     
                        snps_sample_size = sp.sum(nonsyn_snp_filter)
                           
                        if sp.sum(syn_snp_filter)>1:
                            norm_snps = all_norm_snps[syn_snp_filter]
                            positions = all_positions[syn_snp_filter]
                            
                            sample_indices = sorted(sp.random.choice(sp.arange(len(positions)), snps_sample_size, replace=False))
                            norm_snps = norm_snps[sample_indices]
                            positions = positions[sample_indices]
                            M,N = norm_snps.shape
                            
                            ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
                            assert M==len(positions), 'A bug detected.'
                            for i in range(M-1):
                                for j in range(i+1,M):
                                    dist = positions[j] - positions[i]
                                    if dist<max_dist:
                                        ld_dist_dict['syn'][dist]['r2_sum']+=ld_mat[i,j]**2
                                        ld_dist_dict['syn'][dist]['snp_count']+=1.0            
            else:
                #Filtering SNPs with small MAFs
                freqs = g['codon_snp_freqs'][...]
                mafs = sp.minimum(freqs,1-freqs)
                maf_filter = mafs>min_maf
                if sp.sum(maf_filter)>1:
                    is_synonimous_snp = g['is_synonimous_snp'][...]
                    is_nonsynonimous_snp = sp.negative(is_synonimous_snp)
                    syn_snp_filter = is_synonimous_snp*maf_filter
                    nonsyn_snp_filter = is_nonsynonimous_snp*maf_filter
    
                    if sp.sum(syn_snp_filter)>sp.sum(nonsyn_snp_filter):
                        all_norm_snps = g['norm_codon_snps'][...]
                        all_positions = g['codon_snp_positions'][...]
                        norm_snps = all_norm_snps[maf_filter]
                        positions = all_positions[maf_filter]
                        M,N = norm_snps.shape
                        
                        ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
                        assert M==len(positions), 'A bug detected.'
                        for i in range(M-1):
                            for j in range(i+1,M):
                                dist = positions[j] - positions[i]
                                if dist<max_dist:
                                    ld_dist_dict['all'][dist]['r2_sum']+=ld_mat[i,j]**2
                                    ld_dist_dict['all'][dist]['snp_count']+=1.0
        
                        
                        if sp.sum(nonsyn_snp_filter)>1:
                            norm_snps = all_norm_snps[nonsyn_snp_filter]
                            positions = all_positions[nonsyn_snp_filter]
                            M,N = norm_snps.shape
                            
                            ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
                            assert M==len(positions), 'A bug detected.'
                            for i in range(M-1):
                                for j in range(i+1,M):
                                    dist = positions[j] - positions[i]
                                    if dist<max_dist:
                                        ld_dist_dict['nonsyn'][dist]['r2_sum']+=ld_mat[i,j]**2
                                        ld_dist_dict['nonsyn'][dist]['snp_count']+=1.0
                                     
                        snps_sample_size = sp.sum(nonsyn_snp_filter)
                           
                        if sp.sum(syn_snp_filter)>1:
                            norm_snps = all_norm_snps[syn_snp_filter]
                            positions = all_positions[syn_snp_filter]
                            
                            sample_indices = sorted(sp.random.choice(sp.arange(len(positions)), snps_sample_size, replace=False))
                            norm_snps = norm_snps[sample_indices]
                            positions = positions[sample_indices]
                            M,N = norm_snps.shape
                            
                            ld_mat = sp.dot(norm_snps,norm_snps.T)/float(N)
                            assert M==len(positions), 'A bug detected.'
                            for i in range(M-1):
                                for j in range(i+1,M):
                                    dist = positions[j] - positions[i]
                                    if dist<max_dist:
                                        ld_dist_dict['syn'][dist]['r2_sum']+=ld_mat[i,j]**2
                                        ld_dist_dict['syn'][dist]['snp_count']+=1.0
    
    
    
    for plot_type in ld_dist_dict.keys():
        avg_r2s = []
        dist_0_r2s = []
        dist_1_r2s = []
        dist_2_r2s = []
        plot_distances = []
        dist_0s = []
        dist_1s = []
        dist_2s = []
        for dist in distances:
            if ld_dist_dict[plot_type][dist]['snp_count']>0:
                avg_r2 = ld_dist_dict[plot_type][dist]['r2_sum']/float(ld_dist_dict[plot_type][dist]['snp_count'])
                avg_r2s.append(avg_r2)
                plot_distances.append(dist)
                if dist%3==0:
                    dist_0_r2s.append(avg_r2)
                    dist_0s.append(dist)
                elif dist%3==1:
                    dist_1_r2s.append(avg_r2)
                    dist_1s.append(dist)
                elif dist%3==2:
                    dist_2_r2s.append(avg_r2)
                    dist_2s.append(dist)
            
        plot_distances = sp.array(plot_distances)
        dist_0s = sp.array(dist_0s)
        dist_1s = sp.array(dist_1s)
        dist_2s = sp.array(dist_2s)
        avg_r2s = sp.array(avg_r2s)
        dist_0_r2s = sp.array(dist_0_r2s)
        dist_1_r2s = sp.array(dist_1_r2s)
        dist_2_r2s = sp.array(dist_2_r2s)
    
        
        pylab.clf() 
        bins = sp.arange(0,max(plot_distances),bin_size)
        digitize = sp.digitize(plot_distances, bins)    
        xs = []
        ys = []
        for bin_i in range(len(bins)):
            bin_filter = digitize==(bin_i+1)
            if len(plot_distances[bin_filter])>0:
                xs.append(sp.mean(plot_distances[bin_filter]))
                ys.append(sp.mean(avg_r2s[bin_filter]))
        
        pylab.plot(xs, ys, color='k', linestyle='None', marker='.', alpha=0.5)
        pylab.xlabel(r'Pairwise distance ($d$)')
        pylab.ylabel(r'Squared correlation ($r^2$)')
        if filter_pop is not None:
            pylab.savefig('%s/ld_%s_codons_%s.png'%(fig_dir,plot_type,filter_pop))
        else:
            pylab.savefig('%s/ld_%s_codons.png'%(fig_dir,plot_type))
    
        pylab.clf()
        bins = sp.arange(0,max(dist_1s),bin_size)
        digitize = sp.digitize(dist_1s, bins)    
        xs = []
        ys = []
        for bin_i in range(len(bins)):
            bin_filter = digitize==(bin_i+1)
            if len(dist_1s[bin_filter])>0:
                xs.append(sp.mean(dist_1s[bin_filter]))
                ys.append(sp.mean(dist_1_r2s[bin_filter]))
        pylab.plot(xs,ys, linestyle='None', marker='.', color='green', alpha=0.5, label=r'$d$ mod $3 = 1$')
    
        bins = sp.arange(0,max(dist_2s),bin_size)
        digitize = sp.digitize(dist_2s, bins)    
        xs = []
        ys = []
        for bin_i in range(len(bins)):
            bin_filter = digitize==(bin_i+1)
            if len(dist_2s[bin_filter])>0:
                xs.append(sp.mean(dist_2s[bin_filter]))
                ys.append(sp.mean(dist_2_r2s[bin_filter]))    
        pylab.plot(dist_2s,dist_2_r2s, linestyle='None', marker='.', color='red', alpha=0.5, label=r'$d$ mod $3 = 2$')
    
        bins = sp.arange(0,max(dist_0s),bin_size)
        digitize = sp.digitize(dist_0s, bins)    
        xs = []
        ys = []
        for bin_i in range(len(bins)):
            bin_filter = digitize==(bin_i+1)
            if len(dist_0s[bin_filter])>0:
                xs.append(sp.mean(dist_0s[bin_filter]))
                ys.append(sp.mean(dist_0_r2s[bin_filter]))    
        pylab.plot(dist_0s,dist_0_r2s, linestyle='None', marker='.', color='blue', alpha=0.5, label=r'$d$ mod $3 = 0$')
        pylab.xlabel(r'Pairwise distance ($d$)')
        pylab.ylabel(r'Squared correlation ($r^2$)')
        pylab.legend()
        if filter_pop is not None:
            pylab.savefig('%s/part_ld_%s_codons_%s.png'%(fig_dir,plot_type,filter_pop))
        else:
            pylab.savefig(fig_dir+'/part_ld_%s_codons2.png'%(plot_type))

 
def parse_pop_map(file_name = '/project/NChain/faststorage/rhizobium/ld/Rhizobium_soiltypes_nod.txt'):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, origin in izip(t['Seq ID'], t['Genospecies']):
        pop_map[str(strain_id)]=origin
    return pop_map, sp.array(t['Genospecies'])
    

def gen_sfs_plots(snps_hdf5_file = '/project/NChain/faststorage/rhizobium/ld/called_snps.hdf5', 
                 fig_dir = '/project/NChain/faststorage/rhizobium/ld', filter_pop=None):
    pop_map, ct_array = parse_pop_map()
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
                        if pop_map[s]==filter_pop:
                            indiv_filter[s_i]=True
                    except:
                        continue
                if sp.sum(indiv_filter)<20:
                    continue
                codon_snps = g['codon_snps'][...]
                codon_snps = codon_snps[:,indiv_filter]
                t_codon_snps = sp.transpose(codon_snps)
                freqs = sp.mean(t_codon_snps,0)
                
            else:
                freqs = g['codon_snp_freqs'][...]
            mafs = sp.minimum(freqs,1-freqs)
            is_synonimous_snp = g['is_synonimous_snp'][...]
            syn_mafs.extend(mafs[is_synonimous_snp])
            nonsyn_mafs.extend(mafs[sp.negative(is_synonimous_snp)])
            all_mafs.extend(mafs)
            
    
    if filter_pop is not None:
        pylab.clf()
        pylab.hist(all_mafs, bins=50)
        pylab.title('SFS (all binary codon SNPs)')
        pylab.savefig('%s/sfs_all_%s.png'%(fig_dir,filter_pop))
    
        pylab.clf()
        pylab.hist(nonsyn_mafs, bins=50)
        pylab.title('SFS (non-synonimous SNPs)')
        pylab.savefig('%s/sfs_non_syn_%s.png'%(fig_dir,filter_pop))
    
        pylab.clf()
        pylab.hist(syn_mafs, bins=50)
        pylab.title('SFS (synonimous SNPs)')
        pylab.savefig('%s/sfs_syn_%s.png'%(fig_dir,filter_pop))
        
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
        
