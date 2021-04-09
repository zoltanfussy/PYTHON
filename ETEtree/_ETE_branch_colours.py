#make a pallette of colours for a set of locations
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import csv
import os

#################
### FUNCTIONS ###
#################

def find_supported(node, support):
    "Find nodes with significant support"
    matches = []
    for n in node.traverse():
        if n.support >= support and n.is_leaf() == False and n.is_root() == False: 
            matches.append(n)
            n.add_features(name="SUP")
    return matches

def tag_replace(string):
    #print("Lineage lookup began, please use the following list of unrecognized genera to update your fetch_lineages.tsv")
    if "-" in string:
        string = string.replace("-", "_")
    if " " in string:
        string = string.replace(" ", "_")
    tag = string.split("_")[0]
    if tag in {"Candidatus", "uncultured", "mt", "pt", "cyto", "SEED", "SEED1", "SEED2", "WEIRD", "candidate"}:
        tag = string.split("_")[1]
    if tag in taxarepl9:
        genus = taxarepl9.get(tag, "").split(" ")[0]
        string = string.replace(tag, taxarepl9[tag])
    else:
        genus = tag
    if tag in fetch_lineages:
        string = "{}@{}".format(string, fetch_lineages[tag])
    elif genus in fetch_lineages:
        string = "{}@{}".format(string, fetch_lineages[genus])
    else:
        print("No lineage available for: {}".format(tag))
        all_lineages_found = False

    return string
    
#################
### DATA READ ###
#################

filename = "hapto_slow0.9_newick"
prunefile = "hapto_slow0.9_skip.txt"

if os.path.isdir("/Users/morpholino/OwnCloud/"):
    home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
    home = "/Volumes/zoliq data/OwnCloud/"
else:
    print("Please set a homedir")
if True:
    print("setting default directory")
    defdir = "progs/PYTHON/ETEtree/"
    coredata = home + "progs/PYTHON/"
    wd = home + defdir

#load additional data
print("loading data files")  
LocDataFile = csv.reader(open(wd + "4pred-preds.txt"), delimiter="\t", skipinitialspace=True)

taxareplfile = csv.reader(open(coredata + "taxarepl9.tsv"), delimiter="\t", skipinitialspace=False)
taxarepl9 = {row[0]:row[1] for row in taxareplfile}

lineagefile = csv.reader(open(coredata + "fetch_lineages.tsv"), delimiter="\t", skipinitialspace=False)
fetch_lineages = {row[0]:row[1] for row in lineagefile}

all_lineages_found = True
#################
print("loading tree file")  
t = Tree(filename + ".treefile", format=0) #format flexible with support values

#set the ancestors here...
ancestors_d = {"Clade1": ["cyanANABv_WP_011317903.1_plastoquinolTOX", 
                          "grnPRASc_MMETSP0941_Gene.14464-Transcript_5625_Chlorophyta_Prasinococcales"],
               "Clade2": ["crypGONIp_MMETSP0107_Gene.30083-Transcript_20766_Cryptophyta_Cryptomonadales",
                          "Phenylobacterium_sp._RIFCSPHIGHO2_01_FULL_69_31_OHB27812.1"]
              }
#or make them available in a tsv:
with open('_ancestors.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    try:
        ancestors_d.update({r[0]: [r[1], r[2]] for r in reader})
    except IndexError:
        #incomplete line
        print("Incomplete processing of ancestor lineages!")
        print(ancestors_d)

try:
    ancestor = t.get_common_ancestor(ancestors_d[filename])
    #print(ancestor)
    t.set_outgroup(ancestor)
except KeyError:
    print("Root not selected!")
    print(t.get_tree_root())
    quit()
t.ladderize(direction=1)

#select scale 0-1.0 or 0-100 for support values
supportscache = t.get_cached_content(store_attr="support")
supportslist = [x.support for x in supportscache]
if max(supportslist) == 1:
    minsupport = 0.85
else:
    minsupport = 85
find_supported(t, support=minsupport) #find non-terminal nodes with high support

    
##################
###    MAIN    ###
##################
#create a dictionary of taxon data and a list of all localities
locdata = {}
localization = set()

for row in LocDataFile:
    locdata[row[0]] = row[2]
    localization.add(row[2])

#print(location)
list(localization).sort()

#assign colors to localities
N = len(localization)
colors = {"MT": "dodgerblue", "PT": "mediumseagreen", "dual": "darkviolet", 
          "CS": "black", "amb": "black", 
          "SP": "black", "no pred": "black"}
#automatic assignment:
#cmap = get_cmap(N)
#colors = {}
#for i, X in enumerate(localization):
#    colors[X] = cmap[i] #hopefully this is still recognized as a color
#print(colors)


#set different graphic styles:
othersupport, fullsupport, highsupport, leaf = NodeStyle(), NodeStyle(), NodeStyle(), NodeStyle()
fullsupport["hz_line_width"] = 5
fullsupport["size"] = 0
#fullsupport["extra_branch_line_type"] = 2 #the extra branches need to be adjusted still

highsupport["hz_line_width"] = 5
highsupport["hz_line_color"] = "#666666"
highsupport["size"] = 0
#highsupport["extra_branch_line_type"] = 2 #the extra branches need to be adjusted still
othersupport["size"] = 0

#define tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.min_leaf_separation = 1
#ts.show_branch_length = True
ts.show_branch_support = True
ts.scale =  200 # 100 pixels per branch length unit
ts.min_leaf_separation = 0.5
ts.branch_vertical_margin = 0 # 10 pixels between adjacent branches
ts.title.add_face(TextFace(filename, fsize=20), column=0)
circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.scale = 100
circular_style.arc_start = -90 # 0 degrees = 3 o'clock
circular_style.arc_span = 330

pruned = []
#now for something completely different, tree traverse
for node in t.traverse():
    if node.is_leaf():
        pruned.append(node)
        local = locdata.get(node.name, "no pred")
        leafcolor = colors.get(local, "black")
        node.add_features(color=leafcolor)
        node.img_style["size"] = 0
        node.img_style["fgcolor"] = leafcolor
        node.img_style["vt_line_color"] = leafcolor
        node.img_style["hz_line_color"] = leafcolor
        node.img_style["hz_line_width"] = 2
        #node.name = tag_replace(node.name) #make s
        #print(node.name)
        #node.set_style(leaf)
        #print(node.name, local)
    elif node.name == "SUP":
        if node.support == 100:
            node.set_style(fullsupport)
        else:
            node.set_style(highsupport)
    else:
        node.set_style(othersupport)
#to check tree is annotated:
#print(t.get_ascii(attributes=["name", "abundance"], show_internal=True)) 
#color and abundance are new features
#pruned.sort()
print("all terminal branches:", len(pruned))
with open(prunefile) as f:
    skipnodes = f.read().split("\n")
    skipnodes = [x.split("[&!")[0] for x in skipnodes]
    print("To be pruned:\n{}\n{}\n".format(len(skipnodes), "\n".join(skipnodes)))
pruned = [x for x in pruned if x.name not in skipnodes]
print("remain after pruning:", len(pruned))
t.prune(pruned, preserve_branch_length=True)

#ts.legend.add_face(fig, column=0) #does not work
#t.show(tree_style=ts)
t.write(format=0, outfile="{}_pruned.treefile".format(filename))
t.render(filename + "_c.pdf", w=400, units="mm", tree_style=ts)

if all_lineages_found == False:
    print("run fetch_entrez_lineage_tree.py to find missing lineages")

print("Analysis finished")