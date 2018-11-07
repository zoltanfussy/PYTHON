#make a pallette of colours for a set of locations
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import csv
import matplotlib.pyplot as plt
import matplotlib
import math

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
    tag = string.split("_")[0]
    if tag in taxarepl9:
        string = string.replace(tag, taxarepl9[tag])
    return string
    
#################
### DATA READ ###
#################
filename = "Clade1"
#################
t = Tree(filename + ".treefile", format=0) #format flexible with support values
t.ladderize(direction=1)

#set the ancestors here
ancestors_d = {"Clade1": ["cyanANABv_WP_011317903.1_plastoquinolTOX", 
                          "grnPRASc_MMETSP0941_Gene.14464-Transcript_5625_Chlorophyta_Prasinococcales"],
               "Clade1": ["crypGONIp_MMETSP0107_Gene.30083-Transcript_20766_Cryptophyta_Cryptomonadales",
                          "Phenylobacterium_sp._RIFCSPHIGHO2_01_FULL_69_31_OHB27812.1"]
              }
try:
    ancestor = t.get_common_ancestor(ancestors_d[filename])
    t.set_outgroup(ancestor)
except KeyError:
    print("Root not selected!")
    print(t.get_tree_root())
    quit()

#select scale 0-1.0 or 0-100 for support values
supportscache = t.get_cached_content(store_attr="support")
supportslist = [x.support for x in supportscache]
if max(supportslist) == 1:
    minsupport = 0.85
else:
    minsupport = 85
find_supported(t, support=minsupport) #find non-terminal nodes with high support

#load additional data    
LocDataFile = csv.reader(open("4pred-preds.txt"), delimiter="\t", skipinitialspace=True)

taxareplfile = csv.reader(open("taxarepl9.tsv"), delimiter="\t", skipinitialspace=False)
taxarepl9 = {}
for row in taxareplfile:
    taxarepl9[row[0]] = row[1]
    
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
print(colors)


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


#now for something completely different, tree traverse
for node in t.traverse():
    if node.is_leaf():
        local = locdata.get(node.name, "no pred")
        leafcolor = colors.get(local, "black")
        node.add_features(color=leafcolor)
        node.img_style["size"] = 0
        node.img_style["fgcolor"] = leafcolor
        node.img_style["vt_line_color"] = leafcolor
        node.img_style["hz_line_color"] = leafcolor
        node.img_style["hz_line_width"] = 2
        node.name = tag_replace(node.name) #make s
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


#ts.legend.add_face(fig, column=0) #does not work
#t.show(tree_style=ts)
t.render(filename + "_c.png", w=400, units="mm", tree_style=ts)