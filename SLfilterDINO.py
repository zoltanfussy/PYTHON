#skript na extrakci sekvenci s dino SL
print("This script filters an input fasta file based on a presence of a SL sequence anywhere within 50 nt from the beginning or end of a given sequence.")
print("Defining a --SL is optional, otherwise dinoflagellate SL is used.")
print("Note that full length SL 'CCGTAGCCATTTTGGCTCAAG' seems too stringent, so 'TTTGGCTCAAG' is default here.")
print("Suggested batch usage: for i in `ls *.fasta`;do python SLfilter.py -in $i -sl TTTGGCTCAAG;done")

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-in', '--infile', help='fasta to be processed', required=True)
parser.add_argument('-sl', '--SL', help='sequence of the spliced leader', default='TTTGGCTCAAG')

args = parser.parse_args()

inFile = SeqIO.parse(args.infile, 'fasta')
mySL = coding_dna = Seq(args.SL, IUPAC.ambiguous_dna)
revSL = mySL.reverse_complement()
out = open('output-' + args.infile, 'w')

for sequence in inFile:
    name = sequence.name
    seq = sequence.seq
    count = 0
    if mySL in sequence.seq[:50] or revSL in sequence.seq[-50:]:
    	out.write('>{}\n{}\n'.format(name,seq))
    	count += 1

out.close()

print("Filtering done. {} sequences with SL.".format(count))
