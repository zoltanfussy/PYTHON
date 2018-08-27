#function to read fasta and create a list of sequences
print("This script generates a fasta file with shorter fasta headers that needs to be analyzed using signalP and targetP.")
print("The output files of signalP and targetP, along with a phylobayes tree of these sequences are then imported.")
print("The script then outputs the renamed treefile with full taxa names and the prediction added at the end.")
print("Please, enter the prefix (-p) of files to be analyzed...(xxx.fasta, xxx.chain.con.tre)")
print("Optionally, provide a taxa replacement key file in a .tsv format (-a)")

import argparse
import os

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-p', '--prefix', help='FASTA prefix', required=True)
parser.add_argument('-a', '--accessions', help='Accession rename key file', default='predefined')

args = parser.parse_args()

prefix = args.prefix
accessions = args.accessions

print("FILE prefix defined: %s" %(prefix))
if accessions != 'predefined':
	print("Taxa replacement key defined: %s" %(accessions))
else:
	print("Only predefined codes are used.")


def read_fasta(filename):
	infile = open(filename)
	lines = infile.read()
	sequences = lines.split(">")[1:]
	infile.close()
	return sequences

#this part reads and returns signalP predictions
def read_signalP(filename):
	infile = open(filename)
	lines = infile.read()
	predictions = lines.split("\n")[2:]
	infile.close()

	signalPout = {}

	for item in predictions:
		if len(item.split()) == 12:
			#pozor handle je omezeny na 30 znaku, aby to nedelalo bugr v dictionary
			handle = item.split()[0][:30]
			pred = item.split()[9]
			if pred == 'N':
				prediction = 'cytosol'
			elif pred == 'Y':
				prediction = 'signal'
			signalPout[handle] = prediction
	return signalPout

def read_targetP(filename):
	infile = open(filename)
	lines = infile.read()
	predictions = lines.split("\n")[7:-3]
	infile.close()

	targetPout = {}

	for item in predictions:
		if len(item.split()) == 8:
			#pozor handle je osekany na nekolik znaku, aby to nedelalo bugr v dictionary
			handle = item.split()[0]
			pred = item.split()[6]
			if pred == 'C':
				prediction = 'plastid'
			elif pred == 'M':
				prediction = 'mito'
			elif pred == 'S':
				prediction = 'signal'
			else:
				prediction = 'cytosol'
			targetPout[handle] = prediction
	return targetPout

#read predefined taxa codes
taxa = {
	'crypGUILt': 'Guillardia theta', 'rhiAMMOs': 'Ammonia sp.', 'dinGYMNc': 'Gymnodinium catenatum', 'funASPEf': 'Aspergillus fumigatus', 
	'kytORYZs': 'Oryza sativa Japonica', 'amoPAULc': 'Paulinella chromatophora', 'dinLINGp': 'Lingulodinium polyedra', 'oomPHYra': 'Phytophtora ramorum', 
	'hwal': 'Haloquadratum walsbyi', 'alfaNITRh': 'Nitrobacter hamburgensis', 'chryHETak': 'Heterosigma akashiwo', 'mlep': 'Mycobacterium leprae', 
	'chlaLOTam': 'Lotharella amoebiformis', 'strTHRAs': 'Thraustochytrium sp.', 'rhiPLASb': 'Plasmodiophora brassicae', 
	'strNANNO': 'Nannochloropsis gaditana B-31', 'redTIMSo': 'Timspurckia oligopyrenoides', 'chryCHATs': 'Chattonella subsalsa', 
	'rhiROSAs': 'Rosalina sp.', 'tgon': 'Toxoplasma gondii', 'ddis': 'Dictyostelium discoideum', 'rhiSORIs': 'Sorites sp.', 
	'strCAFca': 'Cafeteria Caron', 'mmar': 'Methanococcus maripaludi', 'perkPmari': 'Perkinsus marinus', 'dhan': 'Debaryomyces hansenii', 
	'crei': 'Chlamydomonas reinhardtii', 'hsap': 'Homo sapiens', 'halo': 'Halobacterium sp.', 'gamaYERSp': 'Yersinia pestis', 'sent': 'Salmonella enterica', 
	'firmLISTm': 'Listeria monocytogenes', 'bant': 'Bacillus anthracis', 'osat': 'Oryza sativa Japonica', 'grnBATHp': 'Bathycoccus prasinos', 
	'apiCRYPm': 'Cryptosporidium muris', 'redGALDs': 'Galdieria sulphuraria', 'actiMYCOt': 'Mycobacterium tuberculosis', 'dinHETtr': 'Heterocapsa triquestra', 
	'cele': 'Caenorhabditis elegans', 'cbot': 'Clostridium botulinum', 'rhiPARTg': 'Partenskyella glossopodia', 'dinCRYPc': 'Crypthecodinium cohnii', 
	'cilFAVta': 'Favella taraikaensis', 'deth': 'Dehalococcoides ethenogenes', 'atha': 'Arabidopsis thaliana', 'cilTETRt': 'Tetrahymena thermophila', 
	'cjej': 'Campylobacter jejuni', 'ypes': 'Yersinia pestis', 'archPYROa': 'Pyrobaculum aerophilum', 'redRHSOm': 'Rhodosorus marinus', 
	'redCHONc': 'Chondrus crispus', 'archSULFt': 'Sulfolobus tokodaii', 'cneg': 'Cryptococcus neoformans', 'strECTsi': 'Ectocarpus siliculosus', 
	'cyanSYNEs': 'Synechococcus sp.', 'alfaRICKc': 'Rickettsia conori', 'grnPOLYp': 'Polytomella parva', 'wend': 'Wolbachia endosymbiont', 
	'gamaSHEWb': 'Shewanella baltica', 'oomPYTHu': 'Pythium ultimum', 'chryPTERd': 'Pteridomonas danica', 'aaeo': 'Aquifex aeolicus', 
	'cyanPROCm': 'Prochlorococcus marinus', 'tpal': 'Treponema pallidum', 'redCYANm': 'Cyanidioschyzon merolae', 'haptPRYMp': 'Prymnesium parvum', 
	'betaCUPRn': 'Cupriavidus necator', 'dinCERAf': 'Ceratium fusus', 'cpne': 'Chlamydophila pneumoniae', 'haptPLEUc': 'Pleurochrysis carterae', 
	'strAURAl': 'Aurantiochytrium limacinum', 'glauGLOwi': 'Gloeochaete wittrockiana', 'apiGREGn': 'Gregarina niphandrodes', 
	'dinAMPHm': 'Amphidinium massartii', 'cilEUPfo': 'Euplotes focardii', 'chlaLOTgl': 'Lotharella globosa', 'dinAMPHc': 'Amphidinium carterae', 
	'rhiMINch': 'Minchinia chitonis', 'grnPYRpa': 'Pyramimonas parkeae', 'tbru': 'Trypanosoma brucei', 'vcho': 'Vibrio cholerae', 
	'cyanTHERe': 'Thermosynechococcus elongatus', 'gsul': 'Geobacter sulfurreducens', 'dmel': 'Drosophila melanogaster', 'redPORpu': 'Porphyra purpurea', 
	'grnCHLAr': 'Chlamydomonas reinhardtii', 'cyanCYANp': 'Cyanothece sp.', 'dinOXYma': 'Oxyrrhis marina', 'atum': 'Agrobacterium tumefaciens', 
	'funLACCb': 'Laccaria bicolor', 'excNAEgr': 'Naegleria gruberi', 'tcru': 'Trypanosoma cruzi', 'strTHAps': 'Thalassiosira pseudonana', 
	'mtub': 'Mycobacterium tuberculosis', 'dinALEXc': 'Alexandrium catenella', 'excPERco': 'Percolomonas cosmopolitus', 'dinALEXa': 'Alexandrium andersonii', 
	'cilLITOp': 'Litonotus pictus', 'alfaMAGNm': 'Magnetospirillum magneticum', 'grnCHLOv': 'Chlorella variabilis', 'cpos': 'Coccidioides posadasii', 
	'strAPLAs': 'Aplanochytrium stocchinoi', 'redRHOma': 'Rhodella maculata', 'grnOSTRt': 'Ostreococcus tauri', 'eugEUTn': 'Eutreptiella gymnastica NIES-381', 
	'grnVOLVc': 'Volvox carteri f. nagariensis', 'aniMONOb': 'Monosiga brevicollis', 'chryOCHRO': 'Ochromonas sp.', 'sfle': 'Shigella flexneri', 
	'eugEUTc': 'Eutreptiella gymnastica CCMP1594', 'firmSTAPa': 'Staphyllococcus aureus', 'dinKARmi': 'Karlodinium micrum', 'sman': 'Schistosoma mansoni', 
	'chlaLOToc': 'Lotharella oceanica', 'cilPARAt': 'Paramecium tetraurelia', 'dinDURIb': 'Durinskia baltica', 'apiCHROv': 'Chromera velia', 
	'kytPHYSp': 'Physcomitrella patens', 'grnCOCCs': 'Coccomyxa subellipsoidea', 'cilFABRs': 'Fabrea salina', 'spom': 'Schizosaccharomyces pombe', 
	'ylip': 'Yarrowia lipolytica', 'dinPYROb': 'Pyrodinium bahamense', 'strCAFro': 'Cafeteria roenbergensis', 'amoACANc': 'Acanthamoeba castellanii', 
	'kytARABt': 'Arabidopsis thaliana', 'dinKARb': 'Karenia brevis', 'ggal': 'Gallus gallus', 'betaRALSs': 'Ralstonia solanacearum', 
	'crypRHOsa': 'Rhodomonas salina', 'lmon': 'Listeria monocytogenes', 'apiVOROp': 'Voromonas pontica', 'oana': 'Ornithorhynchus anatinus', 
	'actiCORYd': 'Corynebacter diphteriae', 'betaVERMe': 'Verminephrobacter eiseniae', 'rhiELPHm': 'Elphidium margitaceum', 'ctep': 'Chlorobium tepidum', 
	'alfaAZOSs': 'Azospirillum sp.', 'cilTIARf': 'Tiarina fusus', 'chlaBIGna': 'Bigelowiella natans', 'eugEUGlo': 'Euglena longa', 'ppat': 'Physcomitrella patens', 
	'dinNOCTs': 'Noctiluca scintillans', 'haptEMILh': 'Emiliania huxleyi', 'cyanLYNGp': 'Lyngbya sp.', 'cyanANABv': 'Anabaena variabilis', 
	'calb': 'Candida albicans', 'apiTOXOg': 'Toxoplasma gondii', 'cper': 'Clostridium perfringens', 'tmar': 'Thermotoga maritima', 'wsuc': 'Wolinella succinogenes', 
	'syne': 'Synechococcus sp.', 'pfal': 'Plasmodium falciparum', 'ecol': 'Escherichia coli', 'chlaLOTHs': 'Lotharella sp.', 'betaBURKc': 'Burkholderia cenocepacia', 
	'cilPLATm': 'Platyophrya macrostoma', 'eugEUGgr': 'Euglena gracilis', 'redPORae': 'Porphyridium aerugineum', 'rnor': 'Rattus norvegicus', 
	'bcidFLAVc': 'Flavobacterium columnare', 'mmus': 'Mus musculus', 'glauCYPTg': 'Cyanoptyche gloeocystis', 'mjan': 'Methanocaldococcus jannaschii', 
	'grnPYRam': 'Pyramimonas amylifera', 'ssol': 'Sulfolobus solfataricus', 'cyanNODUs': 'Nodularia spumigena', 'grnDUNte': 'Dunaliella tertiolecta', 
	'actiSTREc': 'Streptomyces coelicolor', 'glauCYApa': 'Cyanophora paradoxa', 'chryVAUCH': 'Vaucheria litorea', 'grnMICRp': 'Micromonas pusilla', 
	'bcidBACTf': 'Bacteroides fragilis', 'apiTHEIe': 'Theileria equi', 'eugTRYPb': 'Trypanosoma brucei', 'aory': 'Aspergillus oryzae', 
	'strPHAtr': 'Phaeodactylum tricornutum', 'dinKRYPf': 'Kryptoperidinium foliaceum', 'amoDICTd': 'Dictyostelium discoideum', 'firmBACIa': 'Bacillus anthracis', 
	'saur': 'Staphylococcus aureus', 'dinHETro': 'Heterocapsa rotundata', 'msed': 'Metallosphaera sedula', 'dinGLENf': 'Glenodinium foliaceum', 
	'tvol': 'Thermoplasma volcanium', 'archTHERv': 'Thermoplasma volcanium', 'apiVITbr': 'Vitrella brassicaformis', 'spne': 'Streptococcus pneumoniae', 
	'bcidPORPg': 'Porphyromonas gingivalis', 'redERYTm': 'Erythrolobus madagascarensis', 'bmaa': 'Brugia malayi', 'bmal': 'Burkholderia mallei', 
	'cilCONDm': 'Condylostoma magnum', 'alfaRHODs': 'Rhodobacter sphaeroides', 'eugLEISm': 'Leishmania major', 'ncra': 'Neurospora crassa', 
	'bsui': 'Brucella suis', 'crypGUIth': 'Guillardia theta', 'mmul': 'Macaca mulatta', 'kytSELAm': 'Selaginella moellendorffii', 'dinSYMBs': 'Symbiodinium sp.', 
	'gamaVIBRc': 'Vibrio cholerae', 'rsol': 'Ralstonia solanacearum', 'apiEIMEt': 'Eimeria tenella', 'strSKEma': 'Skeletonema marinoi', 
	'aful': 'Archaeoglobus fulgidus', 'afum': 'Aspergillus fumigatus', 'cbur': 'Coxiella burnetii', 'archPICRt': 'Picrophilus torridus', 
	'strAUREa': 'Aureococcus anophagefferens', 'drad': 'Deinococcus radiodurans', 'egos': 'Eremothecium gossypii', 'klac': 'Kluyveromyces lactis', 
	'cyanCROCw': 'Crocosphaera watsonii', 'ftul': 'Francisella tularensis', 'rbal': 'Rhodopirellula baltica', 'redMADAe': 'Madagascaria erythrocladoides', 
	'apiNEOSc': 'Neospora caninum', 'amoENTAh': 'Entamoeba histolytica', 'crypCRYpa': 'Cryptomonas paramecium', 'ecab': 'Equus caballus', 
	'funCRYPn': 'Cryptococcus neoformans', 'dinPRORm': 'Prorocentrum minimum', 'cimm': 'Coccidioides immitis', 'scer': 'Saccharomyces cerevisiae', 
	'chryDINOB': 'Dinobryon sp.', 'redERYTa': 'Erythrolobus australicus', 'aniDROSm': 'Drosophila melanogaster', 'bcidPREVr': 'Prevotella ruminicola', 
	'strNANNg': 'Nannochloropsis gaditana', 'dinDINak': 'Dinophysis acuminata'
}

###########################################################################
#main body of the code starts here:                                       #
###########################################################################

sequences = read_fasta(prefix + '.fasta')
outfilesequences = open(prefix + '-out.fasta','w')
sequencesdictionary = {}

#here a list of fasta names is created to backcompare outputs from prediction algorithms
#this will create a dictionary of sequences names (values), where keys are 30-character shortenings of these names
for sequence in sequences:
	thirtychar = sequence.split('\n')[0][:30]
	aminoacids = ''.join(sequence.split('\n')[1:])
	sequencesdictionary[thirtychar] = sequence.split('\n')[0]
	newseq = '>{}\n{}\n'.format(thirtychar, aminoacids)
	outfilesequences.write(newseq)

outfilesequences.close()

#addition of specific taxacodes entered by -a to the taxa vocabulary
if accessions != 'predefined':
	infile = open(accessions)
	accessions = infile.read()
	codes = accessions.split("\n")
	infile.close()
	for code in codes:
		taxa[code.split("\t")[0]] = code.split("\t")[1]
else:
	pass

#SignalP and TargetP are needed to run externally
print ("Please, run the following commands: ")
print ("targetp -P %s-out.fasta > %s-targetp.txt" % (prefix, prefix))
print ("signalp -m %s-mature.fasta %s-out.fasta > %s-signalp.txt" % (prefix, prefix, prefix))
print ("targetp -P %s-mature.fasta > %s-bts.txt" % (prefix, prefix))
raw_input("Waiting for SignalP and TargetP predictions, press Return when targetp.txt, signalp.txt and bts.txt files ready.")


if os.path.isfile("%s-targetp.txt" % (prefix)) == False:
	print("-----\nError: targetP outfile missing\n------\n")
#	os.system('targetp -P %s-out.fasta > %s-targetp.txt' % (prefix, prefix))
elif os.path.isfile("%s-signalp.txt" % (prefix)) == False:
	print("-----\nError: signalP outfile missing\n------\n")
#	os.system('signalp -m %s-mature.fasta %s-out.fasta > %s-signalp.txt' % (prefix, prefix, prefix))
elif os.path.isfile("%s-bts.txt" % (prefix)) == False:
	print("-----\nError: bts outfile missing\n------\n")
#	os.system('targetp -P %s-mature.fasta > %s-bts.txt' % (prefix, prefix))
else:
	print("predictions made, continuing")


#this part returns fasta with changed taxacode to species name - for in-house database
#sequencenamelist = []
#for sequence in sequences:
#	code = sequence.split('_')[0]
#	rest = ''.join(sequence.split('_')[1:])
#	thirtychar = sequence.split('\n')[:30]
#	sequencenamelist.append(sequence.split('\n')[0])
#	if code in taxa:
#		newseq = '>{}_{}'.format(taxa[code], rest)
#	else:
#		pass
#	outfilesequences.write(newseq)

#same thing using BioPython
#from Bio import SeqIO
#
#inFile = SeqIO.parse('input.txt', 'fasta')
#
#with open('Split_fasta.txt', 'a') as result:
#    for sequence in inFile:
#        name = sequence.name
#        sequencenamelist.append(name)
#        taxon = name.split('_')[0]
#        restname = ''.join(name.split('_')[1:])
#        seq = sequence.seq
#        result.write('{}_{}\n{}\n'.format(taxon,restname,seq))

#magic happened
#prediction summaries are written to file:
targetpred = read_targetP(prefix + '-targetp.txt')
signalpred = read_signalP(prefix + '-signalp.txt')
btspred = read_targetP(prefix + '-bts.txt')
finalpredictions = {}
outfilepredictions = open(prefix + '-preds.txt','w')
outfilepredictions.write('gene\tfinalpred\t(TargetP/SignalP/BTS)\n')

for item in sequencesdictionary.keys():
	TP = targetpred.get(item, "no pred")
	SP = signalpred.get(item, "no pred")
	BTS = btspred.get(item, "-")
	if (BTS == 'plastid' or BTS == 'mito'):
		final = 'BTS-plastid'
	elif (TP == 'cytosol' and SP == 'cytosol'):
		final = 'cytosol'
	elif (TP == 'mito' and SP == 'cytosol'):
		final = 'mitochondrion'
	elif (TP == 'plastid' and SP == 'cytosol'):
		final = 'plastid'
	#nezaskodilo by mit podminku SP == 'signal': ptz BTS nepokryva moznost SP+"no pred"
	else:
		final = "revise"
	finalpredictions[item]=final
	outfilepredictions.write('{}\t{}\t({}/{}/{})\n'.format(item, final, TP, SP, BTS))

#taxon renaming in the phylobayes tree
strom = open(prefix + '.chain.con.tre', 'r')
strom_line = strom.readline()

for key in sequencesdictionary.keys():
	code = sequencesdictionary[key].split('_')[0]
	rest = '_'.join(sequencesdictionary[key].split('_')[1:4])
	if key in taxa:
		finalname = '{}_{}'.format(taxa[key], finalpredictions[key])
	elif code in taxa:
		finalname = '{}_{}_{}'.format(taxa[code], rest, finalpredictions[key])
	else:
		finalname = sequencesdictionary[key] + '_' + finalpredictions[key]
	strom_line = strom_line.replace(key, finalname)

with open(prefix + '-fin.tre', 'w') as result:
    result.write(strom_line)

print("Hotovo")
