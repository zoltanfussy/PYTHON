import argparse, textwrap
from Bio import SeqIO
from ete3 import Tree


def adjust_lightness(color, amount=0.5):
	#https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
	import matplotlib.colors as mc
	import colorsys
	try:
		c = mc.cnames[color]
	except:
		c = color
	c = colorsys.rgb_to_hls(*mc.to_rgb(c))
	return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def leaf_nodes(infile, infasta):
	if infasta == "":
		currentseqs = set()
	else:
		currentseqs = get_seqnames(infasta)
		finaldataset = get_seqnames("_dataset.fasta")
		currentseqs = {x.replace("@", "_") for x in currentseqs if x in finaldataset}
		print("{} current seqs to highlight".format(len(currentseqs)))
	seen = set()
		
	lineage_d = {}

	with open(infile, "rt") as f, \
		open("itol_labels.txt", "wt") as result:
		result.write(header_colors)

		for l in f:
			l = l.strip().split(",")
			branchname = l[0].replace("@", "_")
			if l[0].split("_")[-1].startswith("x"):
				clustsize = int(l[0].split("_")[-1].replace("x", ""))
				if clustsize > 100:
					thickness = 2.5
				elif clustsize > 30:
					thickness = 2.0
				elif clustsize > 10:
					thickness = 1.5
				else:
					thickness = 1.0
			elif branchname in currentseqs:
				thickness = 2.0
				seen.add(branchname)
			else:
				thickness = 1.0
			if l[1] not in lineage_d:
				lineage_d[l[1]] = []
			if "CZE" in branchname.upper():
				lineage_d[l[1]].append(branchname)
				result.write("{}\tbranch\t{}\tnormal\t{}\n".format(branchname, color_map_nodes.get(l[1], "#979A9A"), thickness))
			elif branchname[:3] == "NAT":
				lineage_d[l[1]].append(branchname)
				result.write("{}\tbranch\t{}\tnormal\t{}\n".format(branchname, color_map_nodes.get(l[1], "#979A9A"), thickness))
			else:
				lineage_d[l[1]].append(branchname)
	if len(currentseqs - seen) > 0:
		print("{} current seqs not found to highlight!".format(len(currentseqs - seen)))
		print(currentseqs - seen)
	return lineage_d


def get_tree_distances(infile, leaf_list):
	t = Tree(infile)
	root = "NC_045512.2"
	if root in t:
		t.set_outgroup(root)
	else:
		raise ValueError("Could not find root!")

	#remove branches not in tree
	leaf_list = [x.replace("@", "_") for x in leaf_list]
	leaf_list = [x for x in leaf_list if x in t]
	leaf_list = [x for x in leaf_list if x not in ignore]
	if len(leaf_list) == 0:
		return 0, 0
	leaf_ids = [t&x for x in leaf_list]
	distances = [x.get_distance(root, topology_only=True) for x in leaf_ids]
	taxon1 = leaf_list[distances.index(min(distances))]
	taxon2 = leaf_list[distances.index(max(distances))]

	return taxon1, taxon2


def get_seqnames(infasta):
	return {x.name for x in SeqIO.parse(infasta, "fasta")}
	

def clades(lineage_d, infile):
	with open(infile, "rt") as f, \
		open("itol_labels.txt", "at") as result, \
		open("itol_texts.txt", "wt") as result2:
		result2.write(header_texts)

		for clade in color_map_ranges:
			if len(lineage_d.get(clade, [])) < 2: #test if list contains leaf names
				continue
			taxon1, taxon2 = get_tree_distances(infile, lineage_d[clade]) #least/most distant
			if taxon1 == 0:
				print("Dropping {}, no leaves in tree ({})".format(clade, ",".join(lineage_d[clade])))
			#if clade not in color_map_ranges:
			#	print(clade, lineage_d[clade])
			result.write("{}|{}\trange\t{}\t{}\n".format(taxon1, taxon2, color_map_ranges.get(clade, "#DCDCDC"), clade))
			result2.write("{}\t{}\t-1\t{}\tbold\t5\t0\n".format(taxon1, clade, color_map_nodes.get(clade, "#979A9A")))


def main():
	global color_map_nodes
	color_map_nodes = {
					"B.1": "#A6CEE3",
					"P.1": "#A6CEE3",
					"P.2": "#A6CEE3",
					"AY.1": "#9053D1",
					"AY.10": "#9053D1",
					"AY.11": "#9053D1",
					"AY.12": "#9053D1",
					"AY.113": "#924FDA",
					"AY.121": "#D4FB79",
					"AY.122": "#D4FB79",
					"AY.125": "#D4FB79",
					"AY.126": "#D4FB79",
					"AY.127": "#D4FB79",
					"AY.16": "#502E75",
					"AY.17": "#9053D1",
					"AY.18": "#9053D1",
					"AY.19": "#9053D1",
					"AY.2": "#9053D1",
					"AY.20.1": "#D4FB79",
					"AY.20": "#9053D1",
					"AY.21": "#9053D1",
					"AY.23": "#9053D1",
					"AY.26": "#4D4D4D",
					"AY.3": "#9053D1",
					"AY.33": "#502E75",
					"AY.34": "#787878",
					"AY.36": "#787878",
					"AY.38": "#9053D1",
					"AY.4.1": "#33A02C",
					"AY.4.2": "#33A02C",
					"AY.4.3": "#33A02C",
					"AY.4.4": "#33A02C",
					"AY.4.5": "#33A02C",
					"AY.4": "#33A02C",
					"AY.41": "#9053D1",
					"AY.42": "#D4FB79",
					"AY.43": "#B2DF8A",
					"AY.46.2": "#9053D1",
					"AY.46.6": "#9053D1",
					"AY.46": "#9053D1",
					"AY.5.2": "#9053D1",
					"AY.5": "#9053D1",
					"AY.6": "#9053D1",
					"AY.7.1": "#9053D1",
					"AY.7.2": "#9053D1",
					"AY.7": "#9053D1",
					"AY.8": "#9053D1",
					"AY.9": "#FF2F92",
					"AY.9.2": "#FF2F92",
					"AY.98.1": "#D4FB79",
					"AY.98": "#D4FB79",
					"AZ.2": "#338CC8",
					"B.1.1.153": "#B15928",
					"B.1.1.318": "#1F78B4",
					"B.1.1.7": "#F3444B",
					"B.1.1": "#FFFF99",
					"B.1.160": "#FDBF6F",
					"B.1.177": "#FB9A99",
					"B.1.258.3": "#999999",
					"B.1.258": "#999999",
					"B.1.351.2": "#FF7F00",
					"B.1.351": "#FF7F00",
					"B.1.525": "#C65911",
					"B.1.617.2": "#CAB2D6",
					"C.36.3": "#DDD100",
					"C.36": "#ABA200",
					"None": "#BBBBBB",
					"Other": "#999999",
					"BA.1": "#FF5DBC",
					"BA.2": "#C90076"
					}

	global color_map_ranges
	color_map_ranges = {
					#"B.1": "#D9EDF8",
					"B.1.1": "#D9EDF8",
					"B.1.258": "#98D3FA",
					"P.1": "#DCF9C3",
					"P.2": "#CFECCD",
					"AY.113": "#B599D3",
					"AY.121": "#D4FB79",
					"AY.122": "#D4FB79",
					"AY.125": "#D4FB79",
					"AY.126": "#D4FB79",
					"AY.127": "#D4FB79",
					"AY.20.1": "#D4FB79",
					"AY.20": "#B599D3",
					"AY.26": "#B8B8B8",
					"AY.36": "#B8B8B8",
					"AY.4": "#97C492",
					"AY.42": "#E4FBAE",
					"AY.43": "#E4FBAE",
					"AY.46": "#B599D3",
					"AY.68": "#FB739B",
					"AY.7.1": "#F69B68",
					"AY.9.2": "#B6A6C8",
					"AY.98.1": "#D4FB79",
					"AY.98": "#D4FB79",
					"B.1.1.1": "#D3D3D3",
					"B.1.1.7": "#FCA4A5", 
					"B.1.160": "#FDCF93",
					"B.1.177": "#FBC6C5",
					"B.1.221": "#D3D3D3",
					"B.1.318": "#AD7FDE",
					"B.1.351": "#FECE9E",
					"B.1.617.2": "#DACEDF",
					"B.1.617": "#DACEDF",
					"C.36.3": "#FFFF99",
					"C.36": "#FFFF99",
					"BA.1": "#FD8FC6",
					"BA.2": "#FE43A1"
				}

	global header_colors
	header_colors = textwrap.dedent("""
					TREE_COLORS
					SEPARATOR TAB
					DATA
					#First 3 fields define the node, type and color
					#Possible node types are:
					#'range': defines a colored range (colored background for labels/clade)
					#'clade': defines color/style for all branches in a clade
					#'branch': defines color/style for a single branch
					#'label': defines font color/style for the leaf label
					#'label_background': defines the leaf label background color
					#for 'clade' and 'branch':
					# field 4 defines the branch style ('normal' or 'dashed')
					# field 5 defines the branch width scale factor 
					# (eg. with value 0.5, branch width for that clade will be 0.5 the standard width)
					""")

	global header_texts
	header_texts = textwrap.dedent("""
					DATASET_TEXT
					SEPARATOR TAB
					DATASET_LABEL	viral clades
					COLOR	#ff0000
					MARGIN	1
					ALIGN_TO_TREE	1
					DATA
					""")
	parser = argparse.ArgumentParser(description='How to use argparse')
	parser.add_argument('-l', '--lineage', help='Lineage report input', default="lineage_report.csv")
	parser.add_argument('-t', '--tree', help='Treefile for processing', default="")
	parser.add_argument('-c', '--current_seqs', help='FASTA of current sequences to highlight', default="")
	#parser.add_argument('-i', '--itol', help='ITOL tree ID', default="")

	args = parser.parse_args()

	global ignore
	ignore = []

	lineage_d = leaf_nodes(args.lineage, args.current_seqs)
	if args.tree != "":
		clades(lineage_d, args.tree)

	#cgi_start = "http://itol.embl.de/batch_downloader.cgi?tree="
	#cgi_settings = "&format=svg&display_mode=2&inverted=0&normal_rotation=0&bootstrap_display=0&align_labels=1&line_width=1"
	#if args.itol != "":
	#	print("Use the following link to export the tree:")
	#	print("{}{}{}".format(cgi_start, args.itol, cgi_settings))


if __name__ == '__main__':
	main()

