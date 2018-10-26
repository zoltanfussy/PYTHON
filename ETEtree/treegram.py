from gram import TreeGram
#from gram import TreeGram,TreeGramRadial,Gram
#with Gram you can place text, symbols, etc.

sometree = 'myTreeFile.tre'
read(sometree)

t = var.trees[0]
#t = func.readAndPop(sometree) #reads the tree and drops the file
t.node('A').br.name = 'label' #set name to a branch, accessed by name - displayed above it
t.node(1).br.uName = 'X' #set name to a branch, accessed by number - displayed below it
nC = t.node('C')
nC.name = 'whatever you want' #direct renaming of a leaf
#the following needs modules from p4:
from p4 import *
nC.name = r'Ab{\textcolor{blue}{cde}}fgh {\textcolor{blue}{\ding{110}}}' #recognizes string to recolor parts of text and adds a square at the end

tg = TreeGram(t)
#tgr = TreeGramRadial(t) #for radial trees

n = t.node('B')
tg.setBranchULabel(n, 'uLabel') #another way using a function, when tg has been initialized
#node names - specified in newick trees? - can be accessed and modified:
#n.label.color='red!50'
#n.label.textSize='tiny' #big enough
#t.doSmartLabels=True #or 'semi'; not sure if not n.label.doSmartLabels, or even tg.doSmartLabels
#tg.fixTextOverlaps=True

tg.baseName = 'myTree'

#conf files to set some user defaults, such as font (def Helvetica)
#other useful things:
#documentFontSize
#pdfViewer
#or set these in individual gram scripts:
tg.documentFontSize = 5.5
#individually:
#tg.leafLabelSize = 'footnotesize' #LARGE compared to document fontsize
#tg.internalNodeLabelSize = 'tiny'
#tg.branchLabelSize = 'scriptsize'

tg.scale=None #use number to set 1.0 scale length in cm - then you can have the same scale for several trees
tg.setScaleBar(length=0.2, xOffset=0.0, yOffset=-0.6) #but mostly default values are okay
tg.yScale=0.7 #this is the spacing between leaves
#tg.widthToHeight=0.67 #a width to height ratio is automatically calculated, but you can set it manually
tg.showNodeNums=False #or set true if you want to modify individual branches/nodes

#there are styles you can apply to branch / leaf / root / bracket label:
tg.render() #activates styleDict
st = tg.styleDict['leaf']
st.textShape = 'itshape' #italics
if 1:
    # Make a new style, and put it in the
    # styleDict, with a name.
    from gram import GramText
    g = GramText("taxonOfInterest")
    g.textShape = 'itshape'
    g.textSize = 'small'
    g.color = 'white'
    g.draw = 'black'
    g.lineThickness = 'very thick'
    g.fill = 'blue!60'
    g.name = 'TOI'
    g.anchor = 'west'
    tg.styleDict[g.name] = g

    # Apply the style to some of the leaves.
    for nNum in [1,4,5,7]:  # and one internal
        n = t.node(nNum)
        n.label.myStyle = 'TOI'
        n.label.anchor = 'north east' #that probably means to lower left from node

tg.setBrokenBranch(1) #draw a broken // branch, but how to set its length?
tg.setbracket(2,3, text='these group together', leftNode=1) #groups two leaves together, leftNode can specify their common branching node and WHAT? maybe draws a rectangle over the group
# only works with nodeNum, which can however be retrieved by:
#some_nodenumber = t.node('B').nodeNum
tg.setBracket(6, 10, text='another group', leftNode=None, rotated=True) #rotate bracket label 90°, leftNode unspecified
#useful example:
"""
g = tg.setBracket(1, 3, text='Bacteria', leftNode=0)
g.fill = 'blue!15'
g = tg.setBracket(5, 7, text='Archaea', leftNode=4)
g.fill = 'orange!20'
g = tg.setBracket(9, 10, text='Eukaryotes', leftNode=8)
g.fill = 'green!30'
"""
t.bracketsLineUp = False #if you have several overlapping brackets which you probably want to avoid

#if you have several trees, read them all at once, and call them one by one
#this allows placing them in arrays:
t = var.trees[1]
tgB = TreeGram(t.dupe(), scale=1.) #you may want not to alter the original tree
#tgB.gX = 4.5 #set X offset of the second tree positive=right
#tgB.gY = -1. #set Y offset of the second tree positive=up
#tg.grams.append(tgB) #add another tree to the first

#try looping with:
#for i in range(0, len(sometree)):
#	t = var.trees[i]

tg.png() #prints PDF and a PNG file
tg.pdfPage() #prints a tree on a single sheet of paper
