from Bio import SeqIO
import os

files = [x for x in os.listdir(".") if x.startswith("scaffolds_2020") and x.endswith("fasta")]

print(files)

for f in files:
	dataset = f.replace(".fasta", ".split.fasta")
	with open(dataset, "w") as result:
		for seq in SeqIO.parse(f, "fasta"):
			c = 1
			remaining = str(seq.seq)
			if len(remaining) > 7500:
				while len(remaining) > 7500:
					sequence = remaining[:5000]
					seqname = "{}p{}_length{}".format(seq.name.split("_length")[0], c, seq.name.split("_length")[1])
					result.write(">{}\n{}\n".format(seqname, sequence))
					remaining = remaining[5000:]
					c += 1
				seqname = "{}p{}_length{}".format(seq.name.split("_length")[0], c, seq.name.split("_length")[1])
				result.write(">{}\n{}\n".format(seqname, remaining))
			elif len(remaining) > 300:
				result.write(">{}\n{}\n".format(seq.name, remaining))