from Bio import SeqIO
import argparse,random

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-f', '--infile', help='Fasta infile', default='haplotype.aln')
parser.add_argument('-c', '--output_count', help='Count of randomly output seqs', type=int, default=1000)

args = parser.parse_args()

compartments = ["any"]
compartments_d = {x: [] for x in compartments}
lineage_d = {}

seq_d = {x.name: x.seq for x in SeqIO.parse(args.infile, "fasta")}
suffix = args.infile.split(".")[-1]
output = args.infile.replace(suffix, "{}.{}".format(args.output_count, suffix))
outputlin = args.infile.replace(suffix, "{}.csv".format(args.output_count, suffix))

with open("lineage_report_haplo.csv", "rt") as f:
	for l in f:
		data = l.strip().split(",")
		seqid = data[0]
		lineage_d[seqid] = l
		if not seqid.startswith("EPI_ISL"):
			continue
		lineage = data[1]
		if lineage not in compartments:
			#compartments.append(seqid)
			#compartments_d[lineage] = [seqid]
			compartments_d["any"].append(seqid)
		else:
			compartments_d[lineage].append(seqid)

with open(output, "wt") as result:#, open(outputlin, "wt") as result2:
	for c in compartments:
		fullset = compartments_d[c]
		print("{} items: {}".format(c, len(fullset)))
		subset = random.sample(fullset, k=args.output_count) #random.choices() would also choose duplicate sequences!
		for item in subset:
			result.write(">{}\n{}\n".format(item, seq_d.get(item, "")))
			#result2.write(lineage_d[item])

with open("popart_region-" + outputlin, "wt") as result,\
	 open("popart_regions.csv") as f:
	for l in f:
		l = l.strip().split(",")
		if l[0] == "":
			result.write(",".join(l) + "\n")
		elif l[0] in subset:
			result.write("{}\n".format(",".join(l)))


with open("popart_lineage-" + outputlin, "wt") as result,\
	 open("popart_lineages.csv") as f:
	for l in f:
		l = l.strip().split(",")
		if l[0] == "":
			result.write(",".join(l) + "\n")
		elif l[0] in subset:
			result.write("{}\n".format(",".join(l)))

