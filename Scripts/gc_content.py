## GC content of pacbio data

def gc_content(sequence):
	sequence = sequence.upper()
	GC = (sequence.count('G') + sequence.count('C'))/float(len(sequence))
	return (GC, len(sequence))

def parse_fasta(filename):
	file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
	file_separe = file.split('>') #spliting each entry by the > 
	#remove before the loop the extra space:
	file_separe.remove('')
	#print file_separe
	header = []
	sequences = {}
	for entry in file_separe:
		seq = entry.splitlines()
		header = seq[0] #these are the first elements of the list 
		seq = ''.join(seq[1:]) #joining the sequences 
		sequences[header] = seq
	return sequences

contigs = parse_fasta('polished_assembly.fasta')

print ' Contig name      |     GC content     | length'
for contig in contigs.keys():
	print (contig, gc_content(contigs[contig]))
