# parsing the genes related to eps LD block and concatenate those that are core genes
import glob

def parse_fasta(filename):
    
    file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
    file_separe = file.split('>') #spliting each entry by the > 
    file_separe.remove('')
    parse_dict = {}
    header = []
    for entry in file_separe:
        seq = entry.splitlines()
        header = seq[0] #these are the first elements of the list 
        
        genome_id = seq[0].split('|')[1]

        #if pacbio_dict[genome][0] == genome_id:
        #    seq = ''.join(seq[1:]) #joining the sequences 
        parse_dict[header] = seq
    return parse_dict


align_file ='/faststorage/project/NChain/rhizobium/Gene_group_Ontology/Blast_group_genes/' 
gene_parse = []

# Concatenate eps genes

eps_genes = ['group6329', 'group4858', 'group4857', 'group4856', 'group4855', 'group6330', 'group4851', 'group5620', 'group5621', 'group5622', 'group5851', 'group5657', 'group5658', 'group5659', 'group5660', 'group5661', 'group5662', 'group4743', 'group4742']


for gene in eps_genes:
	file_name = align_file + gene + '.fna'
    gene_parse.append(parse_fasta(filename = file_name, genome = genome))
