# LD of all 
# 1. Find a reference genome to sort the genes by location
# 2. Calculate the LD of all versus all
# 3. Find possible LD blocks, besides the one we know (nod genes)
# 4. Find gene functions
# 1. Find a reference genome to sort the genes by location


def kinship_versus_corrected_genes(directory = 'C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/'):

    nod_genes = OrderedDict([(4144, 'nodX'), (4143, 'nodN'), (4142, 'nodM'), (4141, 'nodL'), (4140, 'nodE'), (4139, 'nodF'), (4138, 'nodD'), (4137, 'nodA'), (4136, 'nodC'), (4135, 'nodI'), (4134, 'nodJ'), (4129, 'nifB'), (4128, 'nifA'), (4127, 'fixX'), (4126, 'fixC'), (4125, 'fixB'), (4124, 'fixA'), (4123, 'nifH'), (4122, 'nifD'), (4121, 'nifK'), (4120, 'nifE'), (2448, 'rpoB'), (2140, 'recA')])
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
    
    for index1 in range(0, len(genes)):
    	for index2 in range(index1, len(genes)):

    		gene1 = genes[index1]
    		gene2 = genes[index2]

    		# Genes data
        	name1, snps1, strains_1 = (gene1)
        	name2, snps2, strains_2 = (gene2)

            # This works only if we assume that strains_2 and strains_1 are ordered beforehand.  Are they? They are.
            strain_mask_1 = np.in1d(strains_1, strains_2, assume_unique=True)
            fitered_strains_1 = strains_1[strain_mask_1]
            strain_mask_2 = np.in1d(strains_2, fitered_strains_1, assume_unique=True)
            #fitered_strains_2 = strains_2[strain_mask_2]

			# Construct GRM for a singular gene
        	total_snps_1 = snps[strains_mask_1, :]
        	grm_1 = np.divide(np.dot(snps1, snps1.T), snps1.shape[1])

        	flat_grm_1 = grm_1.flatten()
        	norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean()
        	norm_flat_grm1 = flat_grm_1
        	norm_flat_grm1 = norm_flat_grm1 / sp.sqrt(sp.dot(norm_flat_grm1, norm_flat_grm1))
			
			# Gene 2
			# Construct GRM for a singular gene
        	total_snps_2 = snps[strains_mask_2, :]
        	grm_2 = np.divide(np.dot(snps2, snps2.T), snps2.shape[1])

        	flat_grm_2 = grm_2.flatten()
        	norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean()
        	norm_flat_grm2 = flat_grm_2
        	norm_flat_grm2 = norm_flat_grm2 / sp.sqrt(sp.dot(norm_flat_grm2, norm_flat_grm2))
        
        	corr = pearsonr(norm_flat_grm1, norm_flat_grm2)[0]
        	
        	print corr
        	if corr > 0:
            	r_scores.append(corr)
            	name = gene[0][:-4]
            	original_name.append(name)
            	name = name[5:]

	            # Number of snps    
    	        n_snps.append(snps.shape[1])

        	    # Number of members 
            	n_members.append(snps.shape[0])
    return(r_scores)