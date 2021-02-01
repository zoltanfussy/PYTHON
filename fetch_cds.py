from Bio import Entrez,SeqIO
Entrez.email = 'zoltan.fussy@gmail.com'

def force_taxid(accession):
	print("WARNING: missing taxid in input file taxified.x.out, requesting from NCBI server")
	nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
	#print(nucl.read())
	nucl_record = SeqIO.read(nucl, 'gb')
	print(nucl_record)
	print(nucl_record.features)
	cds = [f for f in nucl_record.features if f.type == "CDS"]
	print(cds.protein_id)
	orgn = nucl_record.annotations['organism']
	
	print("Accession retrieved:", orgn)

	return orgn

orgn = force_taxid("CR932753")