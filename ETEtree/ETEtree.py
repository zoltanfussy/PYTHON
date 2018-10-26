#make a pallette of colours for a set of locations
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import csv
import matplotlib.pyplot as plt
import matplotlib
import math

#################
### FUNCTIONS ###
#################

def get_cmap(n, name='viridis'): #hsv for very divergent data?
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    colormap = plt.cm.get_cmap(name, n)
    rgbcolors = []
    for i in range(colormap.N):
        rgb = colormap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        rgbcolors.append(matplotlib.colors.rgb2hex(rgb))
    return rgbcolors

def find_supported(node, support):
    "Find nodes with significant support"
    matches = []
    for n in node.traverse():
        if n.support > support and n.is_leaf() == False and n.is_root() == False: 
            matches.append(n)
            n.add_features(name="SUP")
    return matches

def nodeshape_size(actual):
    node_size = math.log10(actual) // 1 + 1
    return node_size*3

#################
### DATA READ ###
#################
t = Tree("testtree_sup.tre", format=0) #format flexible with support values
t.ladderize(direction=1)
find_supported(t, support=0.85) #find non-terminal nodes with high support


EnvDataFile = csv.reader(open("leaves_data.txt"), delimiter="\t", skipinitialspace=True)

##################
###    MAIN    ###
##################
#create a dictionary of taxon data and a list of all localities
envdata = {}
location = set()
abundances = set()
for row in EnvDataFile:
    envdata[row[0]] = {"Locality": row[1], "Abundance": row[2]}
    location.add(row[1])
    abundances.add(row[2])
minmax = (min(abundances), max(abundances))
#print(location)
list(location).sort()

#assign colors to localities
N = len(location)
cmap = get_cmap(N)
colors = {}
for i, X in enumerate(location):
    colors[X] = cmap[i] #hopefully this is still recognized as a color
print(colors)
    
for key in envdata:
    locality = envdata[key]["Locality"]
    envdata[key].update({"Color": colors[locality]})

#set different graphic styles:
othersupport, fullsupport, highsupport, leaf = NodeStyle(), NodeStyle(), NodeStyle(), NodeStyle()
fullsupport["hz_line_width"] = 5
fullsupport["size"] = 0
highsupport["hz_line_width"] = 5
highsupport["hz_line_color"] = "#666666"
highsupport["size"] = 0
othersupport["size"] = 0

#define tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.min_leaf_separation = 1
#ts.show_branch_length = True
ts.show_branch_support = True
ts.scale =  200 # 100 pixels per branch length unit
ts.branch_vertical_margin = 8 # 10 pixels between adjacent branches
ts.title.add_face(TextFace("Diatom Tree", fsize=20), column=0)
circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.show_branch_support = True
circular_style.scale = 100
circular_style.arc_start = -90 # 0 degrees = 3 o'clock
circular_style.arc_span = 330


#now for something completely different, tree traverse
for node in t.traverse():
    if node.is_leaf():
        node.add_features(color=envdata[node.name]["Color"], 
                          abundance=envdata[node.name]["Abundance"])
        node.img_style["size"] = nodeshape_size(int(envdata[node.name]["Abundance"]))
        node.img_style["fgcolor"] = envdata[node.name]["Color"]
        #node.set_style(leaf)
        print(node.name, envdata[node.name]["Locality"], envdata[node.name]["Abundance"])
    elif node.name == "SUP":
        if node.support == 1:
            node.set_style(fullsupport)
        else:
            node.set_style(highsupport)
    else:
        node.set_style(othersupport)
#to check tree is annotated:
#print(t.get_ascii(attributes=["name", "abundance"], show_internal=True)) 
#color and abundance are new features


#generate a scale bare with colours
fig=plt.figure(figsize=(1.5,N*1.5))
ax=fig.add_subplot(111)   
plt.axis('scaled')
ax.set_xlim([ 0, N])
ax.set_ylim([-0.5, 0.5])
for i in range(N):
    rect = plt.Rectangle((i, -0.5), 1, 1, facecolor=cmap[i])
    ax.add_artist(rect)
ax.set_yticks([])
ax.set_xticks([])
plt.show()

#ts.legend.add_face(fig, column=0) #does not work
#t.show(tree_style=ts)
t.render("mytree_c.png", w=300, units="mm", tree_style=circular_style)
t.render("mytree.png", w=300, units="mm", tree_style=ts)




