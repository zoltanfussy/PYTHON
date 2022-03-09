import argparse

#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='tRNAscan-SE output to process', required=True)
args = parser.parse_args()
infile = args.infile
outfile = infile.split(".")[0] + ".gff3"

codechange = {"Ala": "A", "Arg": "R", "Asn": "N",
"Asp": "D", "Cys": "C", "Glu": "E", 
"Gln": "Q", "Gly": "G", "His": "H",
"Ile": "I", "Leu": "L", "Lys": "K",
"Met": "M", "Phe": "F", "Pro": "P",
"Ser": "S", "Thr": "T", "Trp": "W",
"Tyr": "Y", "Val": "V", "Sup": "W", "Undet": "X"}

c = 1
with open(infile) as f, open(outfile, "w") as result:
	for line in f:
		data = line.split("\t")
		contig = data[0].strip()
		if contig not in ("Sequence", "Name", "--------"):
			start = data[2].strip()
			stop = data[3].strip()
			tRNA = data[4].strip()
			oneletter = codechange[tRNA]
			score = data[8].strip()
			anticodon = data[5].strip()

			if int(data[6]) > 0:
				if start < stop:
					intron_start = int(data[6]) - 1
					intron_end = int(data[7]) + 1
					result.write("{0}\ttRNAScan-SE\tgene\t{1}\t{2}\t{3}\t+\t.\tID=tRNA-{4}{8}_tRNA;Name=trn{9}({5})_gene\n\
{0}\ttRNAScan-SE\ttRNA\t{1}\t{2}\t{3}\t+\t.\tID=tRNA-{4}{8}_tRNA;Name=tRNA-{4};anticodon={5}\n\
{0}\ttRNAScan-SE\texon\t{1}\t{6}\t{3}\t+\t.\tID=tRNA-{4}{8}_exon;Note=contains predicted Intron\n\
{0}\ttRNAScan-SE\texon\t{7}\t{2}\t{3}\t+\t.\tID=tRNA-{4}{8}_exon;Parent=tRNA-{4}{8}_exon\n"
.format(contig, start, stop, score, tRNA, anticodon.lower(), intron_start, intron_end, c, oneletter))
					c += 1

				else:
					intron_start = int(data[6]) + 1
					intron_end = int(data[7]) - 1
					result.write("{0}\ttRNAScan-SE\tgene\t{1}\t{2}\t{3}\t-\t.\tID=tRNA-{4}{8}_tRNA;Name=trn{9}({5})_gene\n\
{0}\ttRNAScan-SE\ttRNA\t{1}\t{2}\t{3}\t-\t.\tID=tRNA-{4}{8}_tRNA;Name=tRNA-{4};anticodon={5}\n\
{0}\ttRNAScan-SE\texon\t{1}\t{6}\t{3}\t-\t.\tID=tRNA-{4}{8}_exon;Note=contains predicted Intron\n\
{0}\ttRNAScan-SE\texon\t{7}\t{2}\t{3}\t-\t.\tID=tRNA-{4}{8}_exon;Parent=tRNA-{4}{8}_exon\n"
.format(contig, start, stop, score, tRNA, anticodon.lower(), intron_start, intron_end, c, oneletter))
					c += 1
			else:
				if start < stop:
					result.write("{0}\ttRNAScan-SE\tgene\t{1}\t{2}\t{3}\t+\t.\tID=tRNA-{4}{6}_tRNA;Name=trn{7}({5})_gene\n\
{0}\ttRNAScan-SE\ttRNA\t{1}\t{2}\t{3}\t+\t.\tID=tRNA-{4}{6}_tRNA;Name=tRNA-{4};anticodon={5}\n\
{0}\ttRNAScan-SE\texon\t{1}\t{2}\t{3}\t+\t.\tID=tRNA-{4}{6}_exon\n"
.format(contig, start, stop, score, tRNA, anticodon.lower(), c, oneletter))
					c += 1
				else:
					result.write("{0}\ttRNAScan-SE\tgene\t{1}\t{2}\t{3}\t-\t.\tID=tRNA-{4}{6}_tRNA;Name=trn{7}({5})_gene\n\
{0}\ttRNAScan-SE\ttRNA\t{1}\t{2}\t{3}\t-\t.\tID=tRNA-{4}{6}_tRNA;Name=tRNA-{4};anticodon={5}\n\
{0}\ttRNAScan-SE\texon\t{1}\t{2}\t{3}\t-\t.\tID=tRNA-{4}{6}_exon\n"
.format(contig, stop, start, score, tRNA, anticodon.lower(), c, oneletter))
					c += 1

print("Converting done.")


