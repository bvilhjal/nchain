# Identify genes that do not follow the overall population structure, i. e. genes that introgress.
# For this purpose we need to:
# 1. Construct the phylogenetic tree based on all the genes
# 2. Use the correct population structure genes and correlates it to the kinship matrix
import h5py
import numpy as np
import scipy as sp
import pandas as pd
import os
import glob 
import math
from scipy.stats.stats import pearsonr
import pylab as pl
from numpy import linalg
from collections import OrderedDict
from sys import argv

if argv[1] == 'windows':
    snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'
    out_dir='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps_test/'
    meta_data = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'

if argv[1] == 'mac':
    snps_file='/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5'
    out_dir='/Users/PM/Desktop/PHD_incomplete/Methods/Intergenic_LD/corrected_snps_test/'
    meta_data = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'

def parse_pop_map(file_name = meta_data):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, sara_id, origin, country in izip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country}
    return pop_map

def kinship_all_genes(snps_file= snps_file,
                 min_maf=0.10,
                 max_strain_num=198):
    """
    Calculates the kinship. Our expectation is that will be close to identity matrix
    """
    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()
    all_strains = set()
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) == max_strain_num:
            all_strains = set(strains).union(all_strains)
    num_strains = len(all_strains)
    print 'Found %d "distinct" strains' % num_strains

    ordered_strains = sorted(list(all_strains))

    strain_index = pd.Index(ordered_strains)
    K_snps = sp.zeros((num_strains, num_strains))
    counts_mat_snps = sp.zeros((num_strains, num_strains))

    for i, gg in enumerate(gene_groups):
        if i % 100 == 0:
            print 'Working on gene nr. %d' % i
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(strains) < max_strain_num:
            strain_mask = strain_index.get_indexer(strains)

            # Already normalized snps
            snps = data_g['norm_snps'][...]
            freqs = data_g['freqs'][...]
            mafs = sp.minimum(freqs, 1 - freqs)
            maf_mask = mafs > min_maf
            snps = snps[maf_mask]

            if len(snps) == 0:
                continue
            K_snps_slice = K_snps[strain_mask]
            K_snps_slice[:, strain_mask] += sp.dot(snps.T, snps)
            K_snps[strain_mask] = K_snps_slice
            counts_mat_snps_slice = counts_mat_snps[strain_mask]
            counts_mat_snps_slice[:, strain_mask] += len(snps)
            counts_mat_snps[strain_mask] = counts_mat_snps_slice

    K_snps = K_snps / counts_mat_snps  # element-wise division
    print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))

    tuple_index_kinship = (K_snps, ordered_strains)

    headers = list()
    maps = parse_pop_map()
    for i in ordered_strains:
        headers.append( 'SM' +maps[i]['sara_id'])

    np.savetxt('strains_order.csv', headers,  delimiter=",", fmt='%s')
    np.savetxt('kinship_maf_01.csv', K_snps, fmt='%.18e', delimiter=',', header = ','.join(headers), comments="")
    # Saving as a compressed numpy file:
    #np.savez_compressed("{}/{}".format(kinship_matrix), matrix=snps, strains=strains, maf = maf) 

    return(tuple_index_kinship)

def plot_dirty_PCA( figure_fn = 'pca.png', k_figure_fn = 'kinship_heatmap.png', title=None,
                   figure_dir = '/project/NChain/faststorage/rhizobium/ld/figures',strains=None):
    from scipy import linalg
    
    kinship_mat, strains = (kinship_all_genes())

    evals, evecs = linalg.eig(kinship_mat)  #PCA via eigen decomp
    evals[evals<0]=0
    sort_indices = sp.argsort(evals,)
    ordered_evals = evals[sort_indices]
    print ordered_evals[-10:]/sp.sum(ordered_evals)
    pc1,pc2 = evecs[:,sort_indices[-1]],evecs[:,sort_indices[-2]]
    pl.clf()
    

    if strains is not None:    
        ct_marker_map = {'DK':'*','UK':'^', 'F':'o'}
        gs_color_map = {'gsA':'#386CB0','gsB':'#FB8072', 'gsC':'#1B9E77', 'gsE': '#F0027F', 'gsD': '#984EA3'}
        pop_map = parse_pop_map()

        print pop_map
        for i, strain in enumerate(strains):
            print i
            d = pop_map.get(strain,'NA')
            if d=='NA':
                gs = 'NA'
                country = 'NA'
            else:
                gs = d['genospecies']
                country = d['country']
            pl.scatter(pc1[i],pc2[i], marker=ct_marker_map[country], c=gs_color_map[gs], alpha=0.3, s=100, edgecolor='none')
        for gs in gs_color_map:
            pl.scatter([], [], color=gs_color_map[gs], marker = 's', label=gs, s=100, edgecolor='none')
        for country in ct_marker_map:
            if country !='NA':
                pl.scatter([], [], color='k', marker = ct_marker_map[country], label=country, s=100, facecolors='none')

        
        pl.legend(scatterpoints=1)
    
    # Ploting Principal components    
    else:
        pl.plot(pc1,pc2,'k.')
    if title is not None:
        pl.title(title)
    pl.xlabel('PC1')
    pl.xlabel('PC2')
    pl.tight_layout()
    pl.savefig('test',format='pdf')
    pl.clf()
    pl.imshow(kinship_mat, cmap='hot', interpolation='nearest')

    # Ploting the variance explained by each PC 
    tot = sum(evals)
    var_exp = [(i / tot)*100 for i in sorted(evals, reverse=True)]
    np.savetxt('test.out', var_exp, delimiter=',') 
    cum_var_exp = np.cumsum(var_exp)
    print cum_var_exp

    # Ploting the cumulative variance explained
    with pl.style.context('seaborn-whitegrid'):
        pl.figure(figsize=(6, 6))
        pl.bar(range(198), var_exp, alpha= 1, align='center', label='individual explained variance')
        pl.step(range(198), cum_var_exp, where='mid', label='cumulative explained variance')
        pl.ylabel('Explained variance ratio')
        pl.xlabel('Principal components')
        pl.legend(loc='best')
        #plt.tight_layout()
        pl.savefig('variance_explained_PC1234')
        pl.show()

    #pylab.savefig(figure_dir+'/'+k_figure_fn)
#plot_dirty_PCA()

def kinship_pseudo_genes(directory = out_dir,
                        num_strains = 198):
    # Changing directory
    os.chdir(directory)

    genes = []
    for f in glob.glob('*.npz'):
        with np.load(directory + f) as data:
            genes.append((f, data["matrix"], data["strains"]))

    K_snps = sp.zeros((num_strains, num_strains))
    counts_mat_snps = sp.zeros((num_strains, num_strains))

    for gene in genes:
        name, snps, strain_mask = (gene)
        K_snps_slice = K_snps[strain_mask]
        K_snps_slice[:, strain_mask] += sp.dot(snps, snps.T)
        K_snps[strain_mask] = K_snps_slice
        counts_mat_snps_slice = counts_mat_snps[strain_mask]
        counts_mat_snps_slice[:, strain_mask] += len(snps)
        counts_mat_snps[strain_mask] = counts_mat_snps_slice

        if len(strain_mask) == 198:
            final_strain_mask = strain_mask 

    K_snps = K_snps / counts_mat_snps  # element-wise division
    print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))

    ##pl.matshow(K_snps)
    #pl.title('Kinship pseudo genes')
    #pl.colorbar()
    #pl.savefig('kinshi_pseudo_genes.png')
    #pl.show()

    return(final_strain_mask, K_snps)


def simple_tracy_widow(matrix, PCS = 5):
    '''Decompose the covariance matrix of each gene and extract the sum over the first eigenvalues (PCs)'''
    evals, evecs = (linalg.eigh(matrix))

    sorted_evals = sorted(evals)
    # This give us an estimator of the variance explained by the gene
    variance_explained = sorted_evals[0:PCS]

    return(sum(variance_explained))

def nucleotide_diversity_JC(array1, array2):
    ''' Equation 3.8 'Molecular Evolution and phylogenetics (Nei and Kumar): correcting for multiple hits'''    
    N = len(array1)
    numDiffs = 0
    for i in xrange(N):
        if array1[i]!=array2[i]:
            numDiffs += 1
    #if (numDiffs == 0):
    #    distance = 0.0
    #else:
    #    distance = -0.75*math.log(1 - (4.0 / 3.0) * (float(numDiffs) / N)) 
    return float(numDiffs) / N


def nucleotide_diversity(snps_hdf5_file = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/new_snps.HDF5', 
                                 seq_file = '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/snps.hdf5',
                                 geno_species=['gsA', 'gsA'], bin_size=0.2,
                                 gt_hdf5_file= '/Users/PM/Desktop/PHD_incomplete/Bjarnicode/snps.hdf5',
                                 distance_method= 'normal', gene = '3097'):


    h5f = h5py.File(gt_hdf5_file)
    print h5f.keys()
    ag = h5f['alignments']

    phi_rho = list()
    g1 = ag[gene]
    for ind1 in zip(g1['nsequences'], g1['strains']): # sequences are in the first entry
        for ind2 in zip(g1['nsequences'], g1['strains']):
            
            # do not calculate the diagonal
            if ind1[1][:4] != ind2[1][:4]:
                phi_rho.append(nucleotide_diversity_JC(ind1[0], ind2[0]))

    # Stats         
    # Calculating the average pairwise differences
    individual_mean = sum(phi_rho)/len(phi_rho)
    print individual_mean
    return individual_mean

#print nucleotide_diversity(gene = '2939')


def kinship_versus_corrected_genes(directory = out_dir):

    nod_genes = OrderedDict([(4144, 'nodX'), (4143, 'nodN'), (4142, 'nodM'), (4141, 'nodL'), (4140, 'nodE'), (4139, 'nodF'), (4138, 'nodD'), (4137, 'nodA'),
    (4136, 'nodC'), (4135, 'nodI'), (4134, 'nodJ'), (4129, 'nifB'), (4128, 'nifA'), (4127, 'fixX'), (4126, 'fixC'), (4125, 'fixB'), (4124, 'fixA'), 
    (4123, 'nifH'), (4122, 'nifD'), (4121, 'nifK'), (4120, 'nifE'), (2448, 'rpoB'), (2140, 'recA')])
    os.chdir(directory)

	# Upload overall kinship matrix
    k_strains, kinship = (kinship_pseudo_genes())

	#This is the gene SNPs matrix
    genes = []
    for f in glob.glob('*.npz'):
        with np.load(directory + f) as data:
	    	genes.append((f, data["matrix"], data["strains"]))

    r_scores = []
    gene_name = []
    original_name = []
    n_snps = []
    tracy_variance = []
    n_diversity = []
    #locations = gene_locations()
    print len(genes)
    for gene in genes:

        name, snps, strains_1 = (gene)
        strains_mask_1 = np.in1d(strains_1, k_strains, assume_unique = True)
        filtered_strains_1 = strains_1[strains_mask_1]

        strain_mask_2 = np.in1d(k_strains, filtered_strains_1, assume_unique=True)
        filtered_strains_2 = k_strains[strain_mask_2]

		# Construct GRM for a singular gene
        grm_1 = np.divide(np.dot(snps, snps.T), snps.shape[0])

        flat_grm_1 = grm_1.flatten()
        norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
        norm_flat_grm1 = flat_grm_1
        norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))
		
		# Overall kinship (which is the identity matrix)
        grm_2 = kinship[strain_mask_2, :]
        grm_2 = grm_2[:, strain_mask_2]

        flat_grm_2 = grm_2.flatten()
        norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
        norm_flat_grm2 = flat_grm_2
        norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))

        # Mantel test: correlation of flat matrices
        corr = pearsonr(flat_grm_1, flat_grm_2)[0]
        
        if corr > 0:
                    # Simple trace-widow: measure of variance
        #print simple_tracy_widow(norm_flat_grm2)
            tracy_variance.append(simple_tracy_widow(grm_1))
            
            r_scores.append(corr)
            name = gene[0][:-4]
            original_name.append(name)
            name = name[5:]
            print name
            n_diversity.append(nucleotide_diversity( gene = str(name)))

            if int(name) in nod_genes.keys():
                print nod_genes[int(name)]
                gene_name.append(nod_genes[int(name)])
            else:    
                gene_name.append(gene[0][:-4])

            # Number of snps    
            n_snps.append(snps.shape[1])

            # Number of members 

            # Finding the gene location
            #if gene[0][:-4] in locations.keys():
                #print locations[gene[0][:-4]]
                #origin.append(locations[gene[0][:-4]])
            #else:
                #origin.append(0)

    # 1 statistic:
    # suming up the square of the values off the diagonal (to be close to the identity matrix) 

    # 2 statistic: PCA analysis
    # Try doing an eigen decomp, where the statistic of interest would be the sum of say ~5 top eigenvalues.

    #
    # Include the plasmid origin
    # Include the gene functionality
    # Include the counts by genospecies
    # Include non-synonymous/ synonymous ratio (signals of selection, any other parameter?)
    # Count the number of snps in each gene and make a distribution
    # Colour by number of snps and by locations

    os.chdir('/Users/PM/Desktop/PHD_incomplete/nchain')
    np.save('gene_groups_mantel.npy', original_name)
    LD_stats = pd.DataFrame({'r_scores': r_scores,'gene':gene_name, 'snps': n_snps, 'tracy_variance': tracy_variance, 'pi' : n_diversity})
    LD_stats.to_csv('introgressed_gene_stats_test_2.csv', header = True)

kinship_versus_corrected_genes()