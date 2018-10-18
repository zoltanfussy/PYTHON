import os
from Bio import Entrez,SeqIO
Entrez.email = "zoltan.fussy@google.com" 

#pick if you want to work with rename files or fastas:
filetype =  ["fasta", "renamefile"][0]

badchars = ("|@+,:;()'")
allowed = ("fasta", "fas", "fst", "phy", "phylip")
os.chdir("/Users/zoliq/ownCloud/Jankoviny/Tick_transcriptome/fucosyltree/family-specific")
if filetype == "renamefile":
	files = [x for x in os.listdir(".") if x.startswith("rename")]
	for file in files:
		newfile = "new" + file
		print("Processing file " + file)
		with open(file) as f, open(newfile, "w") as out:
			for l in f:
				if l.startswith("gi_"):
					gid = l.split("_")[1]
					handle = Entrez.efetch(db="protein", id=gid, rettype="fasta", retmode="XML")
					record = Entrez.read(handle)
					try:
						accession = record[0]['TSeq_accver']
					except KeyError:
						accession = l.split("\t")[0].replace("gi_{}_".format(gid), "")
					#print(record)
					annot = "".join([x for x in record[0]['TSeq_defline'] if x not in badchars])
					#print(annot)
					out.write("{}\t{}_{}\n".format(l.split("\t")[0], accession, annot))
					#record['TSeq_taxid'] and record['TSeq_orgname'] also possible
				else:
					out.write(l)

elif filetype == "fasta":
	files = [x for x in os.listdir(".") if x.endswith("fasta")]
	#test purposes only: files = ["STT3_pfam02516.fasta"]
	for file in files:
		newfile = "new" + file
		print("Processing file " + file)
		f = SeqIO.parse(file, "fasta")
		with open(newfile, "w") as out:
			for l in f:
				if l.name.startswith("gi|"):
					gid = l.name.split("|")[1]
					handle = Entrez.efetch(db="protein", id=gid, rettype="fasta", retmode="XML")
					record = Entrez.read(handle)
					try:
						accession = record[0]['TSeq_accver']
					except KeyError:
						accession = l.name.replace("gi|{}|".format(gid), "")
					#print(record)
					annot = "".join([x for x in record[0]['TSeq_defline'] if x not in badchars])
					#print(annot)
					out.write(">{} {}\n{}\n".format(accession, annot, l.seq))
					#record['TSeq_taxid'] and record['TSeq_orgname'] also possible
				else:
					description = "".join(x for x in l.description if x not in badchars)
					out.write(">{}\n{}\n".format(description, l.seq))

"""
Swissprot works in a similar way:
>>> from Bio import ExPASy,SwissProt

>>> handle = ExPASy.get_sprot_raw(hitid)
>>> record = SwissProt.read(handle)
>>> dir(record)
['__doc__', '__init__', '__module__', 'accessions', 'annotation_update',
'comments', 'created', 'cross_references', 'data_class', 'description',
'entry_name', 'features', 'gene_name', 'host_organism', 'keywords',
'molecule_type', 'organelle', 'organism', 'organism_classification',
'references', 'seqinfo', 'sequence', 'sequence_length',
'sequence_update', 'taxonomy_id']
"""