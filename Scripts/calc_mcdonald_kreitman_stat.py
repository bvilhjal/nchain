import scipy as sp
import h5py 
import os
import matplotlib
matplotlib.use('Agg')
import pylab
# pylab.rcParams['legend.numpoints'] = 1

import collections
import pandas as pd
import cPickle
import gzip

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


 
def parse_pop_map(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, origin, country in izip(t['Seq ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'genospecies':origin, 'country':country}
    #print pop_map
    return pop_map


def call_good_snps(sequence, ok_snps, snp_positions, codon_syn_map=None, ok_seq_filter=None, seq_num_vars=None):
    from itertools import izip
    codon_syn_map = get_codon_syn_map()
    M,N = ok_snps.shape
    snps = []
    nts = []
    
    codon_snps = []
    codon_snp_positions = []
    codons = []
    aacids = []
    is_synonimous_snp =  []
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
            if not (ok_seq_filter[snp_pos+1] and ok_seq_filter[snp_pos+2]):
                continue
            if not(seq_num_vars[snp_pos+1]==1 and seq_num_vars[snp_pos+2]==1):
                continue
            cdseq12 = sequence[snp_pos+1:snp_pos+3]
            codon0 = nt0+cdseq12
            codon1 = nt1+cdseq12
            
        elif codon_pos==1:
            if not (ok_seq_filter[snp_pos-1] and ok_seq_filter[snp_pos+1]):
                continue
            if not(seq_num_vars[snp_pos-1]==1 and seq_num_vars[snp_pos+1]==1):
                continue
            cdseq0 = sequence[snp_pos-1]
            cdseq2 = sequence[snp_pos+1]
            codon0 = cdseq0+nt0+cdseq2
            codon1 = cdseq0+nt1+cdseq2
        
        elif codon_pos==2:
            if not (ok_seq_filter[snp_pos-1] and ok_seq_filter[snp_pos-2]):
                continue
            if not(seq_num_vars[snp_pos-1]==1 and seq_num_vars[snp_pos-2]==1):
                continue
            cdseq01 = sequence[snp_pos-2:snp_pos]
            codon0 = cdseq01+nt0
            codon1 = cdseq01+nt1
        
        assert codon0!=codon1, 'Codons are identical?'
        
        #This appears to be a unique codon change with a dimorphic SNP.
        codons.append([codon0,codon1])
        #print 'These are the codons', codons
        freq = sp.mean(snp,0)
        #print freq
        
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
                
        codon_snps.append(snp)
        codon_snp_positions.append(snp_pos)
        
    return {'snps':snps, 'nts':nts, 'codon_snps':codon_snps, 'codon_snp_positions':codon_snp_positions, 
            'codons':codons, 'aacids': aacids, 'is_synonimous_snp':is_synonimous_snp,
            'num_syn_sites':tot_num_syn_sites, 'num_non_syn_sites':tot_num_non_syn_sites, }


def calc_mcdonald_kreitman_stat(geno_species=['gsB', 'gsC'], min_num_strains=30, min_num_sub_pol=10,
                                gt_hdf5_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/snps.hdf5',
                                fig_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/figures',
                                out_file = 'mk_stats_gsB_gsC.hdf5'):

    ni_stats = []
    pop = parse_pop_map()
    pop_map = pop.keys()
    ct_array = pop.values()
    codon_syn_map = get_codon_syn_map()
    h5f = h5py.File(gt_hdf5_file)
    ag = h5f['alignments']
    gene_groups = sorted(ag.keys())
    num_parsed_genes = 0
    dn_ds_ratio_dict = {}
    oh5f = h5py.File(out_file,'w')
    for gg in gene_groups:
        g = ag[gg]
        
        #0. Check if there is evidence for CNVs/paralogs?
        seq_ids = g['strains']
        strains_list = sp.array(map(lambda x: x.split('-')[0], seq_ids))

        gs_list = sp.array([pop[strain]['genospecies'] for strain in strains_list])

        gs_filters = [sp.in1d(gs_list,[gs]) for gs in geno_species]
        #print gs_filters
        common_filter = sp.zeros((len(gs_list)),dtype='bool8')
        for i in range(len(geno_species)):
            common_filter += gs_filters[i]

        gs_strains_lists = [strains_list[gs_filter] for gs_filter in gs_filters]
        gs_strains = [ ]
        has_paralogs = False
        for gs_strains_list in gs_strains_lists:
            gs_strains = sp.unique(gs_strains_list)
            has_paralogs = len(gs_strains)<len(gs_strains_list)
            if has_paralogs:
                break
        num_strains = []
        for gs_strains_list in gs_strains_lists:
            num_strains.append(len(gs_strains_list))
        num_strains = sp.array(num_strains)
        #print num_strains
        
        if has_paralogs:
            #pass
            print 'Evidence for paralogs/CNVs'
        elif sp.all(num_strains>min_num_strains):
            gs_strains = gs_strains_lists
            all_gs_strains = strains_list[common_filter]

            gs_list = sp.array([pop[strain]['genospecies'] for strain in all_gs_strains])
            gs_filters = [sp.in1d(gs_list,[gs]) for gs in geno_species]
                        
            #1. Filter rows with indels and missing data
            nt_mat = g['nsequences'][...]
            nt_mat = nt_mat[common_filter]
            
            no_gaps_no_missing = sp.all(nt_mat<5,0)
            nt_mat = sp.transpose(nt_mat)
            if sp.sum(no_gaps_no_missing)>5:
                raw_snps = nt_mat[no_gaps_no_missing]
                
                print 'Working on gene group: %s'%gg
                #First calc within genospcies Ka/Ks
                d = {}
                for i, gs in enumerate(geno_species):
                    gs_filter = gs_filters[i]
                    gs_raw_snps = raw_snps[:,gs_filter]
                    
                    num_vars = sp.apply_along_axis(lambda x: len(sp.unique(x)), 1, nt_mat[:,gs_filter])
                    ok_num_vars = sp.apply_along_axis(lambda x: len(sp.unique(x)), 1, gs_raw_snps)
                    const_seq_filter = ok_num_vars==1
                    good_snp_filter = ok_num_vars==2

                    num_bin_snps = sp.sum(good_snp_filter)
                    if num_bin_snps>5:
                        
                        M,N = nt_mat.shape
                        non_gap_positions = sp.arange(M)[no_gaps_no_missing]
                        
                        #3. Identify good SNPs (dimorphic SNPs)
                        ok_snps = gs_raw_snps[good_snp_filter]
                        snp_positions = non_gap_positions[good_snp_filter]
                        assert len(ok_snps)==len(snp_positions), 'A bug detected!'
                        
                        #4. Call good SNPs                        
                        sequences = (g['sequences'][...])[common_filter] 
                        good_snps_dict = call_good_snps(sequences[0], ok_snps, snp_positions, codon_syn_map = codon_syn_map,
                                                        ok_seq_filter = no_gaps_no_missing, seq_num_vars = num_vars)
                        
#                         codon_snps = good_snps_dict['codon_snps']
                        is_synonimous_snp = good_snps_dict['is_synonimous_snp']
                        num_syn_sites = good_snps_dict['num_syn_sites']
                        num_non_syn_sites = good_snps_dict['num_non_syn_sites']
                                                
#                         norm_codon_snps = sp.transpose(codon_snps)
#                         codon_snp_freqs = sp.mean(norm_codon_snps,0)
                        
                        #Calculate dn/ds ratios
                        num_syn_pol = sp.sum(is_synonimous_snp)
                        num_non_syn_pol = len(is_synonimous_snp)-num_syn_pol
                        if num_syn_pol>0:
                            pn_ps_ratio = (num_non_syn_pol/num_non_syn_sites)/(num_syn_pol/num_syn_sites)
                        else:
                            pn_ps_ratio=-1

                        d[gs]={'pn_ps_ratio':pn_ps_ratio, 'num_syn_pol':num_syn_pol, 'num_non_syn_pol':num_non_syn_pol, 
                               'M':len(nt_mat), 'const_seq_filter':const_seq_filter, 'num_syn_sites':num_syn_sites, 
                               'num_non_syn_sites':num_non_syn_sites}
                    else:
                        d[gs]={'pn_ps_ratio':-1, 'num_syn_pol':0, 'num_non_syn_pol':0, 
                               'M':len(nt_mat), 'const_seq_filter':const_seq_filter,
                               'num_syn_sites':0, 'num_non_syn_sites':0}
                
                
                #Get the constrained seq filter for the two genospecies
                gs1 = geno_species[0]
                gs2 = geno_species[1]                
                const_seq_filter1 = d[gs1]['const_seq_filter']
                const_seq_filter2 = d[gs2]['const_seq_filter']
                constrained_seq_filter = const_seq_filter1 * const_seq_filter2
                
                
                #Filter seq_num_var array to the two genospecies considered
                gs_filter = gs_filters[0]+gs_filters[1]
                num_vars = sp.apply_along_axis(lambda x: len(sp.unique(x)), 1, nt_mat[:,gs_filter])


                constr_seq_len = sp.sum(constrained_seq_filter)
                if constr_seq_len>5:
                    constr_seq = raw_snps[constrained_seq_filter]
                    constr_num_vars = sp.apply_along_axis(lambda x: len(sp.unique(x)), 1, constr_seq)
                    constr_bin_snps_filter = constr_num_vars==2
                    num_const_seq_bin_snps = sp.sum(constr_bin_snps_filter)
                    if num_const_seq_bin_snps>5:
                        gs_specific_snps = constr_seq[constr_bin_snps_filter]
                        
                        #Get positions for constrained SNPs
                        non_gap_positions = sp.arange(len(nt_mat))[no_gaps_no_missing]
                        constrained_positions = non_gap_positions[constrained_seq_filter]
                        constrained_snps_positions = constrained_positions[constr_bin_snps_filter]

                        #4. Call good SNPs                        
                        good_snps_dict = call_good_snps(g['sequences'][0], gs_specific_snps, constrained_snps_positions, codon_syn_map=codon_syn_map,
                                        ok_seq_filter = no_gaps_no_missing, seq_num_vars=num_vars)
                        
                        is_synonimous_snp = good_snps_dict['is_synonimous_snp']
                        num_syn_sites = good_snps_dict['num_syn_sites']
                        num_non_syn_sites = good_snps_dict['num_non_syn_sites']
                                                
#                         norm_codon_snps = sp.transpose(codon_snps)
#                         codon_snp_freqs = sp.mean(norm_codon_snps,0)
                        
                        #Calculate dn/ds ratios
                        num_syn_subt = sp.sum(is_synonimous_snp)
                        num_non_syn_subt = len(is_synonimous_snp)-num_syn_subt
                        if num_syn_subt>0:
                            dn_ds_ratio = (num_non_syn_subt/num_non_syn_sites)/(num_syn_subt/num_syn_sites)
                        else:
                            dn_ds_ratio=-1


                        d['%s_%s'%(gs1,gs2)]={'dn_ds_ratio':dn_ds_ratio, 'num_syn_subt':num_syn_subt, 
                                              'num_non_syn_subt':num_non_syn_subt, 
                                              'constr_seq_len':constr_seq_len, 
                                              'num_const_seq_bin_snps':num_const_seq_bin_snps}                        
                        
                    else:
                        print 'No binary variants were found to be specific to either genospecies within the gene.'
                        d['%s_%s'%(gs1,gs2)]={'dn_ds_ratio':-1, 'num_syn_subt':0, 'num_non_syn_subt':0, 
                                              'constr_seq_len':constr_seq_len, 
                                              'num_const_seq_bin_snps':num_const_seq_bin_snps}
 
                else:
                    print 'No sequence was found to be constrained in both genospecies within the gene.'
                    d['%s_%s'%(gs1,gs2)]={'dn_ds_ratio':-1, 'num_syn_subt':0, 'num_non_syn_subt':0, 
                                            'constr_seq_len':constr_seq_len, 
                                            'num_const_seq_bin_snps':0}

                num_syn_pol = d[gs1]['num_syn_pol']+d[gs2]['num_syn_pol']
                num_non_syn_pol = d[gs1]['num_non_syn_pol']+d[gs2]['num_non_syn_pol']
                num_syn_pol_sites = d[gs1]['num_syn_sites']+d[gs2]['num_syn_sites']
                num_non_syn_pol_sites = d[gs1]['num_non_syn_sites']+d[gs2]['num_non_syn_sites']
                
                if num_syn_pol>0:
                    pn_ps_ratio = (num_non_syn_pol/num_non_syn_pol_sites)/(num_syn_pol/num_syn_pol_sites)
                else:
                    pn_ps_ratio = -1
                    
                num_subt = d['%s_%s'%(gs1,gs2)]['num_syn_subt']+d['%s_%s'%(gs1,gs2)]['num_non_syn_subt']
                num_pol = d[gs1]['num_syn_pol']+d[gs1]['num_non_syn_pol'] + d[gs2]['num_syn_pol']+d[gs2]['num_non_syn_pol']
                #Now calculate the neutrality index (MK statistic)
                if d['%s_%s'%(gs1,gs2)]['dn_ds_ratio']>0 and pn_ps_ratio>=0:
                    ni_stat = float(pn_ps_ratio/float(d['%s_%s'%(gs1,gs2)]['dn_ds_ratio']))
                    if num_subt>min_num_sub_pol and num_pol>min_num_sub_pol:
                        print 'Found NI stat to be %0.3f'%ni_stat
                        ni_stats.append(ni_stat)
                else:
                    ni_stat = -1
                
                mk_alpha = 1-ni_stat
                    
                d['%s_%s'%(gs1,gs2)]['ni_stat']=ni_stat
                d['%s_%s'%(gs1,gs2)]['MK_alpha']=mk_alpha
                d['%s_%s'%(gs1,gs2)]['num_subt']=num_subt
                d['%s_%s'%(gs1,gs2)]['num_pol']=num_pol
                dn_ds_ratio_dict[gg]=d
               
                o_gg = oh5f.create_group(gg)
                o_gg.create_dataset('ni_stat',data=ni_stat)
                o_gg.create_dataset('mk_alpha',data=mk_alpha)
                o_gg.create_dataset('num_subt',data=num_subt)
                o_gg.create_dataset('num_pol',data=num_pol)
                o_gg.create_dataset('num_syn_subt',data=num_syn_subt)
                o_gg.create_dataset('num_non_syn_subt', data = num_non_syn_subt)
                #o_gg.create_dataset('pn', data = (num_non_syn_pol+1/num_non_syn_pol_sites+1))
                #o_gg.create_dataset('ps', data = (num_syn_pol+1/num_syn_pol_sites+1)) 
                o_gg.create_dataset('num_non_syn_pol', data = num_non_syn_pol) # from both groups
                o_gg.create_dataset('num_syn_pol', data = num_syn_pol) # from both groups
                o_gg.create_dataset('num_syn_pol_sites', data = num_syn_pol_sites)
                o_gg.create_dataset('num_non_syn_pol_sites', data = num_non_syn_pol_sites)
                
                o_gg.create_dataset('pn_ps_ratio1',data=d[gs1]['pn_ps_ratio'])
                o_gg.create_dataset('pn_ps_ratio2',data=d[gs1]['pn_ps_ratio'])
                o_gg.create_dataset('pn_ps_ratio',data=pn_ps_ratio)
                o_gg.create_dataset('dn_ds_ratio',data=d['%s_%s'%(gs1,gs2)]['dn_ds_ratio'])

                num_parsed_genes +=1
        else:
            pass
#             print 'Too few strains..'
    print 'Parsed %d'%num_parsed_genes
    oh5f.close()
    
    print 'Number of NI stats: %d'%len(ni_stats)
    ni_stats = sp.array(ni_stats)
    ni_stats[ni_stats<0.005]=0.005
    log_nis = sp.log10(ni_stats)
    pylab.hist(log_nis,bins=100)
    pylab.xlabel(r'$\log(NI)_{10}$ (McDonald-Kreitman Neutrality Index)')
    pylab.savefig(fig_dir+'/MK_stats_%s_%s.png'%(geno_species[0],geno_species[1]))


    return  dn_ds_ratio_dict, ni_stats

calc_mcdonald_kreitman_stat(geno_species=['gsB', 'gsC'], min_num_strains=30, min_num_sub_pol=10,
                                gt_hdf5_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/snps.hdf5',
                                fig_dir = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/figures',
                                out_file = 'mk_stats_gsB_gsC.hdf5')