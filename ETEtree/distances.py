from ete3 import Tree
import csv,os

#################
### DATA READ ###
#################

filename = "concat_slowsites_v1-iq_noshalrt"

if os.path.isdir("/Users/morpholino/OwnCloud/"):
    home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
    home = "/Volumes/zoliq data/OwnCloud/"
else:
    print("Please set a homedir")
if True:
    print("setting default directory")
    defdir = "AndyLab/MMETSP-tree/"
    coredata = home + "progs/PYTHON/"
    wd = home + defdir


print("loading data files")
with open(wd + "subtrees_v4_final/headers_GEMs_ids.txt") as f:
	data = f.read().split("\n")
	tips_d = {x.split("\t")[1]: x.split("\t")[0] for x in data}
	tips = list(tips_d.keys())
	print(tips)
	revtips = tips[::-1]
	print(tips)

print("loading tree file")  
t = Tree(wd + "finaltrees/{}.treefile".format(filename), format=0) #format flexible with support values

ancestors_d = {}
#load ancestors:
with open('_ancestors.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    try:
        ancestors_d.update({r[0]: [r[1], r[2]] for r in reader})
    except IndexError:
        #incomplete line
        pass

try:
    ancestor = t.get_common_ancestor(ancestors_d[filename])
    t.set_outgroup(ancestor)
except KeyError:
    print("Root not selected!")
    print(t.get_tree_root())
    quit()
t.ladderize(direction=1)

#print(t)

with open(wd + "finaltrees/concat_slowsites.matrix", "w") as result,\
open(wd + "finaltrees/concat_slowsites_ids.matrix", "w") as result2:
	result.write(".\t{}\n".format("\t".join(revtips)))
	revtipids = [tips_d[t] for t in revtips]
	result2.write(".\t{}\n".format("\t".join(revtipids)))
	for x in tips:
		result.write(x)
		result2.write(tips_d[x])
		X = t&x
		dists = []
		for y in revtips:
			d = X.get_distance(y)
			dists.append(d)
		strdists = [str(x) for x in dists]
		result.write("\t{}\n".format("\t".join(strdists)))
		result2.write("\t{}\n".format("\t".join(strdists)))
