from bioservices.kegg import KEGG

k = KEGG()
path = k.get("K00855")
kdict = k.parse(path)
print(kdict)
help(kdict)
with open("play.out", "wt") as result:
	result.write("\n".join(kdict.keys()))
