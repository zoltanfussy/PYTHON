import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

print("This script parses Chromera velia targeting predictions into a pandas DataFrame to analyze correlations.")
infile = open("chrove-main.csv")
lines = infile.read().split('\n')
infile.close()

pandasindex = ["KEGG","mtorths","nummtorths", "Met", "SP", "ML-2A", "MFM", "mTP", "ASAFind"]
values_dict = {"KEGG": [],"mtorths": [],"nummtorths": [],"Met": [], "SP": [], "ML2A": [], "MFM": [], "mTP": [], "ASAFind": []}
mtorthologs_d = {}
seqlist = []
n = 0

df = {}
for line in lines:
	line = line.split('\t')
	#columns
	#ID		KEGG	mtorthologs	Met		SP 		ML-2A	MFM		mTP		ASAFind
	#[0]	[1]		[2:7]		[8]		[9]		[10]	[11]	[12]	[13]
	n += 1
	#remove empty items from mt ortholog list and fill the mt_orthologs dict
	seqlist.append(line[0])
	mtorthologs = line[2:7]
	while '' in mtorthologs:
		mtorthologs.remove('') #odstraní nedefinované
	mtorthologs_d[line[0]] = mtorthologs
	try:
		values_dict['KEGG'].append(line[1])
		values_dict['mtorths'].append(list(mtorthologs))
		values_dict['nummtorths'].append(len(list(mtorthologs))/6)
		values_dict['Met'].append(line[8])
		values_dict['SP'].append(float(line[9]))
		values_dict['ML2A'].append(float(line[10]))
		values_dict['MFM'].append(float(line[11]))
		values_dict['mTP'].append(float(line[12]))
		values_dict['ASAFind'].append(line[13])
	except ValueError: 
		print("ValueError: ", line[0])
print(n, "lines processed")
#print(values_dict) #check that the dictionary is okay, e.g. desatinna ciarka sposobuje chybu 
df = pd.DataFrame(values_dict, index = seqlist)
print("Data types:\n{}\n\n\n===========".format(df.dtypes))
#print(df.MFM/df.SP)

print("MitoFates(animal)/SignalP correlation:", np.corrcoef(df.MFM, df.SP, rowvar=False))
print("MultiLoc2(animal)/SignalP correlation:", np.corrcoef(df.ML2A, df.SP))
print("TargetP(mitochondrial)/SignalP correlation:", np.corrcoef(df.mTP, df.SP))
print("#mitochondrial protein orthologs/SignalP correlation:", np.corrcoef(df.nummtorths, df.SP))
print("#mitochondrial protein orthologs/MitoFates(animal) correlation:", np.corrcoef(df.nummtorths, df.MFM))
print("#mitochondrial protein orthologs/MultiLoc2(animal) correlation:", np.corrcoef(df.nummtorths, df.ML2A))
print("#mitochondrial protein orthologs/TargetP correlation:", np.corrcoef(df.nummtorths, df.mTP))
#WHY CORRELATE PRODUCES THREE-DIGIT RESULTS???

#generate plots with matplotlib
x = df.MFM
y = df.SP
"""
nullfmt = NullFormatter()

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(1, figsize=(8, 8))

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
axScatter.scatter(x, y)

# now determine nice limits by hand:
binwidth = 0.25
xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
lim = (int(xymax/binwidth) + 1) * binwidth

axScatter.set_xlim((-lim, lim))
axScatter.set_ylim((-lim, lim))

bins = np.arange(-lim, lim + binwidth, binwidth)
axHistx.hist(x, bins=bins)
axHisty.hist(y, bins=bins, orientation='horizontal')

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())

plt.show()
"""


