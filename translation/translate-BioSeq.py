from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

infile = open('filtered.fasta')
line = infile.read()
infile.close()

seqs = line.split('>')[1:]

out = open('translated.fasta', 'w')

for seq in seqs:
	coding_dna = Seq(''.join(seq.split('\n')[1:]), IUPAC.ambiguous_dna)
	range(0,len(coding_dna)-2,3)
	print len(coding_dna)
	protein_seq = coding_dna.translate()
	print protein_seq
	out.write(">%s\n%s\n" % (seq.split('\n')[0], protein_seq))

out.close()
	
