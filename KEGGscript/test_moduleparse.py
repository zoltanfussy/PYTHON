from bioservices.kegg import KEGG
kegg = KEGG()
pathway = kegg.get("ko01230")
dict_data = kegg.parse(pathway)
print(dict_data)
output = open("modules.txt", "w")
modules_dict = {}
for key in dict_data['MODULE'].keys():
	pathway = kegg.get(key)
	module_data = kegg.parse(pathway)
	#print(module_data)
	orthologs = []
	for ortholog in module_data['ORTHOLOGY'].keys():
		data = [ortholog, module_data['ORTHOLOGY'][ortholog]]
		orthologs.append("_".join(data))
		modules_dict[ortholog] = key
	output.write('{}\t{}\t{}\t{}\n'.format(key, module_data['NAME'], module_data['DEFINITION'], "//".join(orthologs)))

