from bioservices.kegg import KEGG

output = open('eclist.txt', 'w')

kegg = KEGG()
pathway = kegg.get('ath00900')
dict_data = kegg.parse(pathway)
# print(dict_data)

# g = x.get('tbr03440:Tb11.01.0910/aaseq')
# print(g)

# res = x.parse_kgml_pathway("tbr03440")
# print(res['entries'][0])

# for key, value in dict_data['GENE'].items():
# 	print(key, value)

# for gene in dict_data['GENE']:
# 	output.write(gene + '\n')

for value in dict_data['GENE'].values():
	EC = value.split('[EC:')[1]
	EC = EC.split(']')[0]
	EC = EC.replace(' ', '\n')
	output.write(EC + '\n')