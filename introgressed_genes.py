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
from scipy.stats.stats import pearsonr
import pylab as pl
from collections import OrderedDict

def parse_pop_map(file_name = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/scripts/Rhizobium_soiltypes_new.txt'):
    from itertools import izip
    
    pop_map = {}
    t = pd.read_table(file_name)
    t = t.rename(columns=lambda x: x.strip())
    for strain_id, sara_id, origin, country in izip(t['Seq ID'], t['Strain ID'], t['Genospecies'], t['Country']):
        pop_map[str(strain_id)]={'sara_id': sara_id, 'genospecies':origin, 'country':country}
    return pop_map

def kinship_all_genes(snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5',
                 min_maf=0.05,
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

    tuple_index_kinship = (K_snps, strain_mask)

    headers = list()
    maps = parse_pop_map()
    for i in ordered_strains:
        headers.append( 'SM' +maps[i]['sara_id'])

    np.savetxt('strains_order.csv', headers,  delimiter=",", fmt='%s')
    np.savetxt('kinship.csv', K_snps, fmt='%.18e', delimiter=',', header = str(headers))
    # Saving as a compressed numpy file:
    #np.savez_compressed("{}/{}".format(kinship_matrix), matrix=snps, strains=strains, maf = maf) 

    return(tuple_index_kinship)
#kinship_all_genes()


def gene_locations( ):

    # Make a dictionary that contains a gene and it respective location, based on the gene scaffold analysis
    # 1. Open the csv files containing the genes and their scaffolds origins 
    # 2. Make a dictionary where each gene will have a specific origin
    # 3. Combine this origin with the mantel results: nod genes versus all ( Mantel_test_nod_all_core)
    df1 = pd.read_csv('C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Gene_locations_0_1000.csv', sep=',', header = 0, index_col = 0)
    df2 = pd.read_csv('C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Gene_locations_1001_2000.csv', sep=',', header = 0, index_col = 0)
    df3 = pd.read_csv('C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Gene_locations_2001_3000.csv', sep=',', header = 0, index_col = 0)
    df4 = pd.read_csv('C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Gene_locations_3001_4000.csv', sep=',', header = 0, index_col = 0)
    df5 = pd.read_csv('C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain/Gene_locations_4001_4307.csv', sep=',', header = 0, index_col = 0)

    frames = [df1, df2, df3, df4, df5]

    df = pd.concat(frames, axis = 1)

    # Making the location dictionary: keys are genes and values are plasmid origin: 
    locations = {}

    for column in df:
        # Taking the mode of each gene location:
        all_locations = df[column].tolist()
        origin = max(set(all_locations), key=all_locations.count)
        #print origin

        if column not in locations:
            locations[column] = origin
        else:
            locations[column].append(origin)
    return locations


def kinship_pseudo_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/',
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

        # Here snps are N x M dimensions
        K_snps_slice = K_snps[strain_mask]
        K_snps_slice[:, strain_mask] += np.dot(snps, snps.T)
        K_snps[strain_mask] = K_snps_slice
        counts_mat_snps_slice = counts_mat_snps[strain_mask]
        counts_mat_snps_slice[:, strain_mask] += snps.shape[1] # Markers present in each gene
        counts_mat_snps[strain_mask] = counts_mat_snps_slice

        if len(strain_mask) == 198:
            final_strain_mask = strain_mask 

    K_snps = K_snps / counts_mat_snps  # element-wise division
    print 'The mean of the GRM diagonal is %f' % np.mean(np.diag(K_snps))

    pl.matshow(K_snps)
    pl.title('Kinship pseudo genes')
    pl.colorbar()
    pl.savefig('kinship_pseudo_genes.pdf')
    pl.show()

    tuple_index_kinship = (K_snps, final_strain_mask)

    return(tuple_index_kinship)

#print kinship_pseudo_genes()

def simple_tracy_widow(matrix, PCS = 5):
    '''Decompose the covariance matrix of each gene and extract the sum over the first eigenvalues (PCs)'''
    evals, evecs = (linalg.eigh(matrix))

    # This give us an estimator of the variance explained by the gene
    variance_explained = evals[PCS]

    return(sum(variance_explained))


def kinship_versus_corrected_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/'):

    nod_genes = OrderedDict([(4144, 'nodX'), (4143, 'nodN'), (4142, 'nodM'), (4141, 'nodL'), (4140, 'nodE'), (4139, 'nodF'), (4138, 'nodD'), (4137, 'nodA'),
    (4136, 'nodC'), (4135, 'nodI'), (4134, 'nodJ'), (4129, 'nifB'), (4128, 'nifA'), (4127, 'fixX'), (4126, 'fixC'), (4125, 'fixB'), (4124, 'fixA'), 
    (4123, 'nifH'), (4122, 'nifD'), (4121, 'nifK'), (4120, 'nifE'), (2448, 'rpoB'), (2140, 'recA')])
    os.chdir(directory)

	# Upload overall kinship matrix
    kinship, k_strains = (kinship_pseudo_genes())

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
    locations = gene_locations()
    print len(genes)
    for gene in genes:

        name, snps, strains_1 = (gene)
        strains_mask_1 = np.in1d(strains_1, k_strains, assume_unique = True)
        filtered_strains_1 = strains_1[strains_mask_1]

        strain_mask_2 = np.in1d(k_strains, filtered_strains_1, assume_unique=True)
        filtered_strains_2 = k_strains[strain_mask_2]

		# Construct GRM for a singular gene
        total_snps_1 = snps[strains_mask_1, :]
        grm_1 = np.divide(np.dot(snps, snps.T), snps.shape[1])
        #pl.matshow(grm_1)
        #pl.show()

        flat_grm_1 = grm_1.flatten()
        #norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
        norm_flat_grm1 = flat_grm_1
        norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))
		
		# Overall kinship
        grm_2 = kinship[strain_mask_2, :]
        grm_2 = grm_2[:, strain_mask_2]

        flat_grm_2 = grm_2.flatten()
        #norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
        norm_flat_grm2 = flat_grm_2
        norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))

        # Mantel test: correlation of flat matrices
        corr = pearsonr(norm_flat_grm1, norm_flat_grm2)[0]

        # Simple trace-widow: measure of variance
        print simple_tracy_widow(norm_flat_grm2)
        print simple_tracy_widow(norm_flat_grm1)
        
        if corr > 0:
            r_scores.append(corr)
            name = gene[0][:-4]
            original_name.append(name)
            name = name[5:]

            if int(name) in nod_genes.keys():
                print nod_genes[int(name)]
                gene_name.append(nod_genes[int(name)])
            else:    
                gene_name.append(gene[0][:-4])

            # Number of snps    
            n_snps.append(snps.shape[1])

            # Number of members 
            n_members.append(snps.shape[0])

            # Finding the gene location
            if gene[0][:-4] in locations.keys():
                #print locations[gene[0][:-4]]
                origin.append(locations[gene[0][:-4]])
            else:
                origin.append(0)

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

    os.chdir('C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain')
    np.save('gene_groups_mantel.npy', original_name)
    LD_stats = pd.DataFrame({'r_scores': r_scores,'gene':gene_name, 'snps': n_snps, 'members': n_members, 'origin': origin})
    LD_stats.to_csv('introgressed_gene_stats_test.csv', header = True)

kinship_versus_corrected_genes()