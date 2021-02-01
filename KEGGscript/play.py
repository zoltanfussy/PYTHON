from bioservices.kegg import KEGG

k = KEGG()
path = k.get("map01100")
kdict = k.parse(path)
print(kdict)
