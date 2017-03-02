
def simple_intergenic_ld_core(max_strain_num=198,
                            maf=0.1,
                            n=10,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Core versus core"""

    t0 = time.time()
    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()

    core_genes = []
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(set(strains)) == max_strain_num:
            core_genes.append(gg)
    core_genes = sorted(core_genes)

    r_scores = []
    gene_names = []
    gene_grm_dict = {}

    # Subsetting core genes
    core_genes = core_genes[0:n]

    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(core_genes), len(core_genes)))
    cor_matrix = pd.DataFrame(cor_matrix, index=core_genes, columns=core_genes)

    for i, gg1 in enumerate(core_genes):
        data_g1 = h5f[gg1]
        total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns
        total_snps_1 = minor_allele_filter(total_snps_1, 0.1)

        ''' 3. Calculate the Kinship matrix for each gene '''
        grm = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])

        flat_grm = grm.flatten()
        norm_flat_grm1 = flat_grm - flat_grm.mean() / sp.sqrt(sp.dot(flat_grm, flat_grm))
        gene_grm_dict[gg1] = {'grm':grm , 'norm_flat_grm':norm_flat_grm1}

        for j, gg2 in enumerate(core_genes):
            if i > j:

                norm_flat_grm2 = gene_grm_dict[gg2]['norm_flat_grm']
                covariance = sp.dot(norm_flat_grm1, norm_flat_grm2)
                var1 = np.sum(abs(norm_flat_grm1 - norm_flat_grm1.mean()) ** 2)
                var2 = np.sum(abs(norm_flat_grm2 - norm_flat_grm2.mean()) ** 2)
                r = covariance / sp.sqrt(var1 * var2)

                cor_matrix[gg2][gg1] = r

                # Checking the values with a scipy built function:
                # r_bel = pearsonr(norm_flat_grm1, norm_flat_grm2)
                # print round(r, 5) == round(r_bel[0], 5)

    cor_matrix.to_csv('Mantel_test_all_all.csv', header=True)
    correlation_plot(cor_matrix)
    t1 = time.time()

    total = t1 - t0
    print 'total amount of time consumed is %f' % total

# simple_intergenic_ld_core()

def simple_intergenic_ld_nod_genes(max_strain_num=198,
                            maf=0.2,
                           snps_file='C:/Users/MariaIzabel/Desktop/MASTER/PHD/Bjarnicode/new_snps.HDF5'):
    """Gives a specific list of genes (nod genes) and calculate LD of these genes with all"""

    nod_genes = nod_genes_f()
    # Decoding the nod gene names
    nod_list = []
    for i in nod_genes.keys():
        nod_list.append(str(i).decode("utf-8"))
    print nod_genes.values()

    h5f = h5py.File(snps_file)
    gene_groups = h5f.keys()

    print h5f[nod_list[0]]
    core_genes = []
    for gg in gene_groups:
        data_g = h5f[gg]
        strains = data_g['strains'][...]
        if len(set(strains)) == max_strain_num:
            core_genes.append(gg)
    core_genes = sorted(core_genes)
    gene_grm_dict = {}

    # Making the correlation matrix to be updated
    cor_matrix = np.zeros((len(nod_list), len(core_genes[0:10])))
    cor_matrix = pd.DataFrame(cor_matrix, index=nod_genes.values(), columns=core_genes[0:100])

    for i, gg1 in enumerate(nod_list):
        try:
            strains_1 = h5f[gg1]['strains'][...]
        except KeyError:
            print 'The nod gene %s is not in our subset of the data' % nod_genes[int(gg1)]
            continue

        for j, gg2 in enumerate(core_genes[0:10]):

            strains_2 = h5f[gg2]['strains'][...]

            set_1, set_2 = set(strains_1), set(strains_2)
            intersec = list(set_1 & set_2)

            strain_mask_2 = []
            strain_mask_1 = []

            for i in intersec:
                strain_mask_2.append(np.unique(strains_2).tolist().index(i))
                strain_mask_1.append(np.unique(strains_1).tolist().index(i))

            strain_mask_2 = sorted(strain_mask_2)
            strain_mask_1 = sorted(strain_mask_1)

            if gg1 not in gene_grm_dict:
                data_g1 = h5f[gg1]
                total_snps_1 = data_g1['snps'][...].T  # strains in the rows, snps in the columns

                # Calculating GRM
                total_snps_1 = minor_allele_filter(total_snps_1, maf)
                grm_1 = np.divide(np.dot(total_snps_1, total_snps_1.T), total_snps_1.shape[1])
                gene_grm_dict[str(gg1)] = {'grm':grm_1}

            if gg2 not in gene_grm_dict:
                data_g2 = h5f[gg2]
                total_snps_2 = data_g2['snps'][...].T  # strains in the rows, snps in the columns

                # Calculating GRM
                total_snps_2 = minor_allele_filter(total_snps_2, maf)
                grm_2 = np.divide(np.dot(total_snps_2, total_snps_2.T), total_snps_2.shape[1])
                gene_grm_dict[str(gg2)] = {'grm':grm_2}

            # Calculating correlation and covariance based on the common subset of strains
            grm_1 = gene_grm_dict[str(gg1)]['grm']
            sub_grm_1 = grm_1[strain_mask_1, strain_mask_1]
            flat_grm_1 = sub_grm_1.flatten()
            norm_flat_grm1 = flat_grm_1 - flat_grm_1.mean() / sp.sqrt(sp.dot(flat_grm_1, flat_grm_1))

            grm_2 = gene_grm_dict[str(gg2)]['grm']
            sub_grm_2 = grm_2[strain_mask_2, strain_mask_2]
            flat_grm_2 = sub_grm_2.flatten()
            norm_flat_grm2 = flat_grm_2 - flat_grm_2.mean() / sp.sqrt(sp.dot(flat_grm_2, flat_grm_2))

            # Built in function, it returns correlation coefficient and the p-value for testing non-correlation
            r = pearsonr(norm_flat_grm1, norm_flat_grm2)
            cor_matrix[gg2][nod_genes[int(gg1)]] += r[0]

    correlation_plot(cor_matrix)
    cor_matrix.to_csv('Mantel_test_nod_all.csv', header=True)
    return cor_matrix

# simple_intergenic_ld_nod_genes()
