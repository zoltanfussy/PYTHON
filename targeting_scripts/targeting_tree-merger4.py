#function to read fasta and create a list of sequences
print("This script accepts prediction outputs, along with a phylobayes tree.")
print("The script then outputs the renamed treefile with full taxa names and the prediction added at the end.")
print("Please, enter the prefix (-p) of prediction files to be analyzed...(xxx.txt)")
print("Treefiles are read automatically from the current folder, or can be entered with (-t).")
print("Optionally, provide a taxa replacement key file in a .tsv format (-a)\n###############################################")

import argparse
import os
import sys
from Bio import SeqIO

def query_yes_no(question, default="yes"):
	"""Ask a yes/no question via raw_input() and return their answer.

	"question" is a string that is presented to the user.
	"default" is the presumed answer if the user just hits <Enter>.
		It must be "yes" (the default), "no" or None (meaning
		an answer is required of the user).

	The "answer" return value is True for "yes" or False for "no".
	"""
	valid = {"yes": True, "y": True, "ye": True,
			 "no": False, "n": False}
	if default is None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise ValueError("invalid default answer: '{}'" .format(default))

	while True:
		sys.stdout.write(question + prompt)
		if sys.platform in ["darwin", "win32"]:
			choice = input().lower()
		elif sys.platform.startswith("linux"):
			choice = raw_input().lower()
		else:
			print("unrecognized OS, check how to use raw input")
			choice = raw_input().lower()

		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "
							 "(or 'y' or 'n').\n")

#### Change to workdir ####
###########################

#homedir = "/Users/zoliq/ownCloud/"
homedir = "/Volumes/zoliq data/ownCloud/"
#wd = homedir + "genomes/phatr/phatr mitoglyco/huge alignments/PASTA alignments nonconverging/goods trimmed/final trees/"
wd = homedir + "genomes/euglena longa/trees/MTOX/RESULT/"
datadir = homedir + "progs/PYTHON/targeting_script/"
os.chdir(wd)

#### Collect Input ####
#######################

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-p', '--prefix', help='Prediction files prefix', default='4pred')
parser.add_argument('-t', '--treefile', help='Treefile', default='none')
parser.add_argument('-a', '--accessions', help='Accession rename key file', default='none')
parser.add_argument('-d', '--directory', help='Working directory', default='.')

args = parser.parse_args()

prefix = args.prefix
treefile = args.treefile
accessions = args.accessions
os.chdir(args.directory)
#ONLY FOR TEST PURPOSES:
#accessions = "leaf_renaming.txt"

print("FILE prefix defined: %s" %(prefix))
if accessions != 'none':
	print("Taxa replacement key defined: %s" %(accessions))
else:
	print("Only predefined codes are used.")

#### Functions ####
###################

def second_largest(numbers):
    count = 0
    m1 = m2 = float('-inf')
    for x in numbers:
        count += 1
        if x > m2:
            if x >= m1:
                m1, m2 = x, m1            
            else:
                m2 = x
    return m2 if count >= 2 else None


#### Open and parse predictions ####
####################################

preds_d = {} #predictions dictionary
leavesfrompreds = set()
#preds_d['accession'] = {'hectar':'ND', 'signalp': 'ND', 'asafind': 'ND', 'targetp': 'ND', 'ML2ANIM': 'ND', 'ML2PLANT': 'ND'}


asafind = open(prefix + "-asafind.txt").read().split('\n')
possiblepredsasafind = {'No idea': 'N/A', 'Yes, Low': 'PLASTID, LOW', 'NA': 'NON-PLASTID', 'No': 'NON-PLASTID', 'Yes, High': 'PLASTID, HIGH'}
for item in asafind:
	#Identifier	SignalP	ASAfind cleavage position	ASAfind/SignalP cleavage site offset	ASAfind 20aa transit score	ASAfind Prediction
	#[0]		[1]		[2]							[3]										[4]							[5]	
	
	if not item.startswith('Identifier') and len(item) != 0:
		item = item.split('\t')
		name = item[0].split("@")[0]
		leavesfrompreds.add(name)
		pred = item[5]
		pred = possiblepredsasafind[pred]
		preds_d[name] = {'asafind': pred}


#HECTAR is not so great
hectar = open(prefix + "-hectar.txt").read().split('\n')
possiblepredshectar = {'mitochondrion': 'MITOCHONDRION', 'no signal peptide or anchor': 'OTHER', 'other localisation': 'OTHER', 'chloroplast': 'PLASTID', 'signal anchor': 'SIGNAL', 'signal peptide': 'SIGNAL'}
for item in hectar:
	#protein id	predicted targeting category	signal peptide score	type II signal anchor score	chloroplast score	mitochondrion score	other score
	#[0]		[1]								[2]						[3]							[4]					[5]					[6]
	#BUT!
	#protein id	step	error message	signal peptide score	type II signal anchor score
	#[0]		[1]		[2]				[3]						[4]

	if not item.startswith('protein id') and len(item) != 0:
		item = item.split('\t')
		item = list(filter(None, item)) #remove list() function in python2.x
		name = item[0].split("@")[0]
		pred = item[1]

		if pred not in ("Failed(U/X)", "error message"):
			#STRAMENOPILES PREDICTIONS:
			#protein id	predicted targeting category	signal peptide score	type II signal anchor score		chloroplast score	mitochondrion score		other score
			#[0]		[1]								[2]						[3]								[4]					[5]						[6]
			
			#OPISTHOKONTS PREDICTIONS:
			#protein id	predicted targeting category	signal peptide score	type II signal anchor score		mitochondrion score		other score
			#[0]		[1]								[2]						[3]								[4]						[5]					

			#UNSPECIFIED EUKARYOTES PREDICTIONS:
			#protein id	predicted targeting category	signal peptide score	type II signal anchor score		other score
			#[0]		[1]								[2]						[3]								[4]
			
			#ERROR:
			#protein id	step	error message	signal peptide score	type II signal anchor score
			#[0]				[1]				[2]						[3]
			pred = possiblepredshectar[pred]
			if len(item) == 7:
				#HECTARPREDICTION "stramenopiles"
				item = [k.replace('-', '0') for k in item]
				hectarpreds_l = [float(item[2]), float(item[4]), float(item[5]), float(item[6])]
				hectarpreds_d = {float(item[2]): 'SP', float(item[4]): 'cTP', float(item[5]): 'mTP', float(item[6]): 'other'}
				pred = ("{}_({}:{} > {}:{})".format(pred, hectarpreds_d[max(hectarpreds_l)], max(hectarpreds_l), hectarpreds_d[second_largest(hectarpreds_l)], second_largest(hectarpreds_l)))
			elif len(item) == 6:
				#HECTARPREDICTION "animals+fungi"
				item = [k.replace('-', '0') for k in item]
				hectarpreds_l = [float(item[2]), float(item[4]), float(item[5])]
				hectarpreds_d = {float(item[2]): 'SP', float(item[4]): 'mTP', float(item[5]): 'other'}
				pred = ("{}_({}:{} > {}:{})".format(pred, hectarpreds_d[max(hectarpreds_l)], max(hectarpreds_l), hectarpreds_d[second_largest(hectarpreds_l)], second_largest(hectarpreds_l)))
			elif len(item) == 6:
				#HECTARPREDICTION "eukaryotes"
				item = [k.replace('-', '0') for k in item]
				hectarpreds_l = [float(item[2]), float(item[4])]
				hectarpreds_d = {float(item[2]): 'SP', float(item[4]): 'other'}
				pred = ("{}_({}:{} > {}:{})".format(pred, hectarpreds_d[max(hectarpreds_l)], max(hectarpreds_l), hectarpreds_d[second_largest(hectarpreds_l)], second_largest(hectarpreds_l)))
			curdict = {'hectar': pred}
			try:
				preds_d[name].update(curdict)
			except KeyError:
				print("Hectar: bad key " + name)
				preds_d[name] = {'hectar': pred}
				leavesfrompreds.add(name)
		else:
			pred = "undefined amino acids U/X"
			try:
				preds_d[name].update(curdict)
			except KeyError:
				print("Hectar: bad key " + name)
				preds_d[name] = {'hectar': pred}
				leavesfrompreds.add(name)


ML2ANIM = open(prefix + "-ML2animal.txt").read().split('\n')
possiblepredsml2 = {"chloroplast": "PLASTID", 'cytoplasmic': "CYTOSOL", "ER": "ER", "extracellular": "EXTRACELLULAR", 
"Golgi apparatus": "GOLGI", "mitochondrial": "MITOCHONDRION", "nuclear": "NUCLEUS", "lysosomal": "LYSOSOME", 
"peroxisomal": "PEROXISOME", "plasma membrane": "PLASMA MEMBRANE", "secretory pathway": "SECRETORY", "vacuolar": "VACUOLE"}
if ML2ANIM[0].startswith("MultiLoc2"):
	ML2ANIM = ML2ANIM[5:]
	for item in ML2ANIM:
		#protein	decreasing localization predictions from: cytoplasmic, ER, extracellular, Golgi apparatus, mitochondrial, nuclear, lysosomal, peroxisomal, plasma membrane
		if len(item) > 0:
			item = item.split('\t')
			name = item[0].split("@")[0]
			Loc = possiblepredsml2[item[1].split(':')[0]]
			pred = ("{}_({} > {})".format(Loc,item[1], item[2]))
			try:
				preds_d[name].update({"ML2ANIMAL": pred})
			except KeyError:
				print("MultiLoc2-animal: bad key " + name)
				preds_d[name] = {'ML2ANIMAL': pred}
				leavesfrompreds.add(name)


ML2PLANT = open(prefix + "-ML2plant.txt").read().split('\n')
if ML2PLANT[0].startswith("MultiLoc2"):
	ML2PLANT = ML2PLANT[5:]
	for item in ML2PLANT:
		#protein	decreasing localization predictions from: chloroplast, cytoplasmic, ER, extracellular, Golgi apparatus, mitochondrial, nuclear, lysosomal, peroxisomal, plasma membrane, vacuolar
		if not item.startswith('#') and len(item) != 0:
			item = item.split('\t')
			name = item[0].split("@")[0]
			Loc = possiblepredsml2[item[1].split(':')[0]]
			pred = ("{}_({} > {})".format(Loc,item[1], item[2]))
			try:
				preds_d[name].update({"ML2PLANT": pred})
			except KeyError:
				print("MultiLoc2-plant: bad key " + name)
				preds_d[name] = {'ML2PLANT': pred}
				leavesfrompreds.add(name)


signalp = open(prefix + "-signalp.txt").read().split('\n') #sensitive version 4.1
for item in signalp:
	# SignalP-3.0 euk predictions  
	# name 	Cmax 	pos ? 	Ymax	pos ? 	Smax 	pos ? 	Smean 	?	D 	? 	name 	!	Cmax	pos ? 	Sprob	?
	#  0	1		2	3	4		5	6	7		8	9	10		11	12	13	14		15	16		17	18	19		20
	# SignalP-4.1 euk predictions  
	# names Cmax	pos 	Ymax	pos 	Smax	pos 	Smean 	D 	?	Dmaxcut 	Networks-used 
	#  0	1		2		3		4		5		6		7		8	9	10			11

	if not item.startswith('#') and len(item) != 0:
		item = item.split()
		name = item[0].split("@")[0]
		if len(item) == 21:
			#SIGNALPVERSION 3.0
			pred = item[13]
			pred = pred.replace("N", "OTHER")
			pred = pred.replace("Y", "SIGNAL")
		if len(item) == 12:
			#SIGNALPVERSION 4.x
			pred = item[9]
			pred = pred.replace("N", "OTHER")
			pred = pred.replace("Y", "SIGNAL")
		try:
			preds_d[name].update({'signalp': pred})
		except KeyError:
			print(item)
			print("SignalP: bad key " + name)
			preds_d[name] = {'signalp': pred}
			leavesfrompreds.add(name)


targetp = open(prefix + "-targetp.txt").read().split('\n') #PLANT + NONPLANT
possiblepredstargetp = {"M": "MITOCHONDRION", "C": "PLASTID", "S": "SIGNAL", "_": "OTHER"}
for item in targetp:
	#Name       Len     mTP     SP      other   Loc     RC
	#[0]		[1]		[2]		[3]		[4]		[5]		[6]
	#BUT PLANT!
	#Name       Len     cTP     mTP     SP      other   Loc     RC
	#[0]		[1]		[2]		[3]		[4]		[5]		[6]		[7]
	#PLANT incl Cleavage Site Prediction:
	#Name       Len     cTP     mTP     SP      other   Loc     RC 		TPlen
	#[0]		[1]		[2]		[3]		[4]		[5]		[6]		[7]		[8]

	if len(item.split()) > 1 and item.split()[1].isnumeric(): #takes only lines with prediction
		item = item.split()
		name = item[0].split("@")[0]
		if len(item) == 7:
			targetpreds_l = [float(item[2]), float(item[3]), float(item[4])]
			targetpreds_d = {float(item[2]): 'mTP', float(item[3]): 'SP', float(item[4]): 'other'}
			Loc = possiblepredstargetp[item[5]]
			pred = ("{}_({}:{} > {}:{})".format(Loc, targetpreds_d[max(targetpreds_l)], max(targetpreds_l), targetpreds_d[second_largest(targetpreds_l)], second_largest(targetpreds_l)))
		elif len(item) == 8 or len(item) == 9:
			targetpreds_l = [float(item[2]), float(item[3]), float(item[4]), float(item[5])]
			targetpreds_d = {float(item[2]): 'cTP', float(item[3]): 'mTP', float(item[4]): 'SP', float(item[5]): 'other'}
			Loc = possiblepredstargetp[item[6]]
			pred = ("{}_({}:{} > {}:{})".format(Loc, targetpreds_d[max(targetpreds_l)], max(targetpreds_l), targetpreds_d[second_largest(targetpreds_l)], second_largest(targetpreds_l)))
		else:
			pred = "nd"
		try:
			preds_d[name].update({'targetp': pred})
		except KeyError:
			print("TargetP: bad key " + name)
			preds_d[name] = {'targetp': pred}
			leavesfrompreds.add(name)
	else:
		pass

print("preds_dictionary collected")


#### Read own codes ####
########################

#read predefined taxa codes
taxarepl9 = {"actiCORYd": "Corynebacter diphteriae", "actiMYCOt": "Mycobacterium tuberculosis", 
"actiSTREc": "Streptomyces coelicolor", "alfaAZOSs": "Azospirillum sp. B506", 
"alfaMAGNm": "Magnetospirillum magneticum", "alfaNITRh": "Nitrobacter hamburgensis", 
"alfaRHODs": "Rhodobacter sphaeroides", "alfaRICKc": "Rickettsia conori", 
"amoACANc": "Acanthamoeba castellanii str neff", "amoDICTd": "Dictyostelium discoideum", 
"amoENTAh": "Entamoeba histolytica", "amoFILAn": "Filamoeba nolandi", "amoPOLYp": "Polysphondylium pallidum pn500", 
"amoSEXAN": "Sexangularia sp. CB-2014", "amoSTYGA": "Stygamoeba regulata", "amoVANEL": "Vannella", 
"apiBABEb": "Babesia bovis", "apiCHROv": "Chromera velia", "apiCRYPm": "Cryptosporidium muris", 
"apiEIMEt": "Eimeria tenella", "apiGREGn": "Gregarina niphandrodes", "apiNEOSc": "Neospora caninum", 
"apiTHEIe": "Theileria equi", "apiTOXOg": "Toxoplasma gondii", "apiVITbr": "Vitrella brassicaformis", 
"apiVORp": "Voromonas pontica", "archPICRt": "Picrophilus torridus", "archPYROa": "Pyrobaculum aerophilum", 
"archSULFt": "Sulfolobus tokodaii", "archTHERv": "Thermoplasma volcanium", "bcidBACTf": "Bacteroides fragilis", 
"bcidFLAVc": "Flavobacterium columnare", "bcidPORPg": "Porphyromonas gingivalis", 
"bcidPREVr": "Prevotella ruminicola", "betaBURKc": "Burkholderia cenocepacia", "betaCUPRn": "Cupriavidus necator", 
"betaRALSs": "Ralstonia solanacearum", "betaVERMe": "Verminephrobacter eiseniae", 
"chlaAMORC": "Amorphochlora amoebiformis", "chlaBIGna": "Bigelowiella natans", "chlaCHLOR": "Chlorarachnion reptans", 
"chlaLOTam": "Lotharella amoebiformis", "chlaLOTgl": "Lotharella globosa", "chlaLOTgZ": "Lotharella globosa", 
"chlaLOTHs": "Lotharella sp. CCMP622", "chlaLOToc": "Lotharella oceanica", "chlaPARTg": "Partenskyella glossopodia", 
"chryCHATs": "Chattonella subsalsa CCMP2191", "chryDINOB": "Dinobryon sp UTEXLB2267", 
"chryHETak": "Heterosigma akashiwo CCMP3107", "chryOCHRO": "Ochromonas sp. CCMP1393", 
"chryVAUCH": "Vaucheria litorea CCMP2940", "cilCLIMv": "Climacostomum virens", "cilMESOD": "Mesodinium pulex", 
"cilLITOp": "Litonotus pictus",
"cilPARAt": "Paramecium tetraurelia", "cilPLATm": "Platyophrya macrostoma", "cilPROTa": "Protocruzia adherens", 
"cilPSEUp": "Pseudocohnilembus persalinus", "cilSTENc": "Stentor coeruleus", "cilTETRt": "Tetrahymena thermophila", 
"crypCHROO": "Chroomonas cf. mesostigmatica", "crypCRYpa": "Cryptomonas paramecium", 
"crypCRYPp": "Cryptomonas paramecium", "crypGONIp": "Goniomonas pacifica", "crypGUILt": "Guillardia theta", 
"crypGUIth": "Guillardia theta", "crypPALPb": "Palpitomonas bilix", "crypRHODl": "Rhodomonas lens", 
"cyanANABv": "Anabaena variabilis", "cyanCROCw": "Crocosphaera watsonii", "cyanCYANp": "Cyanothece sp. PCC 7425", 
"cyanLYNGp": "Lyngbya sp. PCC 8106", "cyanNODUs": "Nodularia spumigena", "cyanPROCm": "Prochlorococcus marinus", 
"cyanSYNEs": "Synechococcus sp. PCC 7335", "cyanTHERe": "Thermosynechococcus elongatus", 
"dinALEXa": "Alexandrium andersonii CCMP2222", "dinALEXc": "Alexandrium catenella OF101", 
"dinALEXt": "Alexandrium tamarense", "dinAMPHc": "Amphidinium carterae", "dinCERAf": "Ceratium fusus PA161109", 
"dinCRYPc": "Crypthecodinium cohnii WH Provasly-Seligo", "dinCRYPZ": "Crypthecodinium cohnii", 
"dinDURIb": "Durinskia baltica CSIRO CS-38", "dinDURIZ": "Durinskia baltica", "dinDINac": "Dinophysis acuminata",
"dinGLENf": "Glenodinium foliaceum CCAP 1116-3", "dinGYMNc": "Gymnodinium catenatum", 
"dinHETro": "Heterocapsa rotundata SCCAP K-0483", "dinHETtr": "Heterocapsa triquestra CCMP 448", 
"dinKARb": "Karenia brevis SP1", "dinKARL": "Karlodinium veneficum", "dinKARmi": "Karlodinium micrum CCMP2283", 
"dinKARZ": "Karenia brevis", "dinKRYPf": "Kryptoperidinium foliaceum", "dinLINGp": "Lingulodinium polyedra CCMP 1738", 
"dinNOCTs": "Noctiluca scintillans", "dinOXYma": "Oxyrrhis marina", "dinPERKm": "Perkinsus marinus", 
"dinPRORm": "Prorocentrum minimum", "dinPYROb": "Pyrodinium bahamense", "dinSCRIP": "Scrippsiella hangoei", 
"dinSYMBs": "Symbiodinium sp.", "eugANGOd": "Angomonas deanei", "eugBODOs": "Bodo saltans", 
"eugEUGgr": "Euglena gracilis", "eugEUGlo": "Euglena longa", "eugEUTc": "Eutreptiella gymnastica CCMP1594", 
"eugEUTn": "Eutreptiella gymnastica NIES-381", "eugLEISm": "Leishmania major", "eugNEOBd": "Neobodo designis", 
"eugPHYTO": "Phytomonas sp isolate em1", "eugTRYPb": "Trypanosoma brucei", "excNAEgr": "Naegleria gruberi", 
"excPERco": "Percolomonas cosmopolitus ATCC50343", "firmBACIa": "Bacillus anthracis", 
"firmLISTm": "Listeria monocytogenes", "firmSTAPa": "Staphyllococcus aureus", "funASPEf": "Aspergillus fumigatus", 
"funCRYPn": "Cryptococcus neoformans", "funDEBAh": "Debaryomyces hansenii cbs767", "funLACCb": "Laccaria bicolor", 
"funNEURc": "Neurospora crassa", "funPUCCg": "Puccinia graminis", "funSCHIp": "Schizosaccharomyces pombe", 
"gamaSHEWb": "Shewanella baltica", "gamaVIBRc": "Vibrio cholerae", "gamaYERSp": "Yersinia pestis", 
"glauGLOEw": "Gloeochaete wittrockiana", "glauCYApa": "Cyanophora paradoxa", 
"glauCYPTg": "Cyanoptyche gloeocystis SAG4.97", "glauGLOwi": "Gloeochaete wittrockiana SAG46.84", 
"grnBATHp": "Bathycoccus prasinos", "grnCHLAl": "Chlamydomonas leiostraca", "grnHELIs": "Helicosporidium sp.",
"grnDUNAt": "Dunaliella tertiolecta", "grnMICpu": "Micromonas pusilla CCMP1545", "grnMICRZ": "Micromonas pusilla", 
"grnNEPHp": "Nephroselmis pyriformis", "grnOSTRm": "Ostreococcus mediterraneus", "grnPICCL": "Picochlorum", 
"grnPICCY": "Picocystis salinarum", "grnPOLYp": "Polytomella parva", "grnPOLYZ": "Polytomella parva", 
"grnPRASc": "Prasinococcus capsulatus", "grnPYRam": "Pyramimonas amylifera CCMP720", 
"grnPYRpa": "Pyramimonas parkeae CCMP726", "grnPYRpZ": "Pyramimonas parkeae", "grnTETRa": "Tetraselmis astigmatica", 
"grnTETRs": "Tetraselmis striata", "grnVOLVc": "Volvox carteri f. nagariensis", "hapCALCl": "Calcidiscus leptoporus", 
"hapEXANg": "Exanthemachrysis gayraliae", "hapIMANT": "Imantonia", "hapISOCHg": "Isochrysis galbana", 
"hapPHEOC": "Phaeocystis", "hapPLEUc": "Pleurochrysis carterae", "hapPRYMp": "Prymnesium parvum", 
"hapSCYPH": "Scyphosphaera apsteinii", "haptEMILh": "Emiliania huxleyi", "haptEMIZ": "Emiliania huxleyi", 
"haptPLEUc": "Pleurochrysis carterae", "haptPRYMp": "Prymnesium parvum Texoma1", 
"hetPERCO": "Percolomonas cosmopolitus", "kytAMBOt": "Amborella trichopoda", "kytARABt": "Arabidopsis thaliana", 
"kytGLYCm": "Glycine max", "kytHORDv": "Hordeum vulgare", "kytORYZs": "Oryza sativa Japonica", 
"kytPHYSp": "Physcomitrella patens", "kytSELAm": "Selaginella moellendorffii", "kytSOlyc": "Solanum lycopersicum", 
"kytZEAma": "Zea mays", "metANOLc": "Anolis carolinensis", "metCAENe": "Caenorhabditis elegans", 
"metDAPHp": "Daphnia pulex", "metHELOr": "Helobdella robusta", "metLOTTg": "Lottia gigantea", 
"metNEMAv": "Nematostella vectensis", "metSCHIm": "Schistosoma mansoni", "metSTROp": "Strongylocentrotus purpuratus", 
"metTRICa": "Trichoplax adhaerens", "opiACANT": "Acanthoeca", "opiMONbr": "Monosiga brevicollis", 
"opiMONOb": "Monosiga brevicollis", "opiSALPr": "Salpingoeca rosetta", "redCHONc": "Chondrus crispus", 
"redCOMPS": "Compsopogon caeruleus", "redCYANm": "Cyanidioschyzon merolae", "redERYTa": "Erythrolobus australicus", 
"redERYTa": "Erythrolobus australicus", "redERYTm": "Erythrolobus madagascarensis", 
"redERYTm": "Erythrolobus madagascarensis", "redGALDs": "Galdieria sulphuraria", 
"redMADAe": "Madagascaria erythrocladioides", "redMADAe": "Madagascaria erythrocladioides", 
"redPORae": "Porphyridium aerugineum SAG 1380-2", "redPORae": "Porphyridium aerugineum", 
"redPORpu": "Porphyra purpurea", "redRHOma": "Rhodella maculata CCMP736", "redRHOSO": "Rhodosorus marinus", 
"redTIMSP": "Timspurckia oligopyrenoides", "rhiAMMOs": "Ammonia sp.", "rhiAMMOZ": "Ammonia sp.", 
"rhiELPHm": "Elphidium margitaceum", "rhiMINch": "Minchinia chitonis", "rhiPARTg": "Partenskyella glossopodia RCC365", 
"rhiPAULc": "Paulinella chromatophora", "rhiPLASb": "Plasmodiophora brassicae", "rhiROSAL": "Rosalina", 
"rhiSORIT": "Sorites sp.", "strALBUl": "Albugo laibachii", "strAPHAi": "Aphanomyces invadans", 
"strAPLAN": "Aplanochytrium", "strAURAl": "Aurantiochytrium limacinum", "strAURAZ": "Aurantiochytrium limacinum", 
"strAUREN": "Aureococcus anophagefferens", "strAUREE": "Aureococcus anophagefferens", 
"strAUREa": "Aureococcus anophagefferens", "strAUREl": "Aureoumbra lagunensis", 
"strBICOS": "Bicosoecida sp. CB-2014", "strBLASh": "Blastocystis hominis", "strBOLID": "Bolidomonas", 
"strBOLIp": "Bolidomonas pacifica", "strCAFro": "Cafeteria roenbergensis", "strCAFsp": "Cafeteria", 
"strCAFca": "Cafeteria str Caron", "strCHATs": "Chattonella subsalsa", "strCHRRH": "Chrysoreinhardia", 
"strCYLIN": "Cylindrotheca closterium", "strDETON": "Detonula confervacea", "strDICTY": "Dictyocha speculum", 
"strDINOB": "Dinobryon", "strDITYL": "Ditylum brightwellii", "strECTsi": "Ectocarpus siliculosus", 
"strEXTUs": "Extubocellulus spinifer", "strFRAGk": "Fragilariopsis kerguelensis", "strHETak": "Heterosigma akashiwo", 
"strMALLO": "Mallomonas", "strNANNg": "Nannochloropsis gaditana B31", "strNITZs": "Nitzschia", 
"strOCHRO": "Ochromonas", "strODONa": "Odontella aurita", "strPARAb": "Paraphysomonas bandaiensis", 
"strPELAG": "Pelagomonas calceolata", "strPHAtr": "Phaeodactylum tricornutum v3", "strPHEOM": "Phaeomonas parva", 
"strPHYpa": "Phytophthora parasitica", "strPHYra": "Phytophtora ramorum", "strPINGU": "Pinguiococcus pyrenoidosus", 
"strPTERd": "Pteridomonas danica", "strPYTHu": "Pythium ultimum var. sporangiiferum BR650", 
"strPYTHv": "Pythium vexans", "strRHIZC": "Rhizochromulina marina", "strSAPRd": "Saprolegnia diclina", 
"strSKEco": "Skeletonema costatum", "strSTAUR": "Staurosira", "strSYNCR": "Synchroma pusillum", 
"strTHAps": "Thalassiosira pseudonana", "strTHNEM": "Thalassionema frauenfeldii", 
"strTHNEn": "Thalassionema nitzschioides", "strTHRAU": "Thraustochytrium sp.", "strTHTRX": "Thalassiothrix antarctica", 
"strVAUCl": "Vaucheria litorea"}


#EXPORT TAXACODES TO FILE
# with open(datadir + "taxarepl9.tsv", "w") as taxaout:
# 	for k in taxarepl9.keys():
# 		taxaout.write("{}\t{}\n".format(k, taxarepl9[k]))
# quit()


#addition of specific taxacodes entered by -a to the taxa vocabulary
if accessions != 'none':
	codes = open(accessions).read().split("\n")
	#codes = open("leaf_renaming.txt").read().split("\n")
	for code in codes:
		code = code.split("\t")	
		if len(code) == 2:
			taxarepl9[code[0]] = code[1]

#genus-to-high taxon conversion file
genus2hightaxonfile = open(datadir + "high_taxon_assignment.txt").read().split("\n")
genus2hightaxon_d = {}
for line in genus2hightaxonfile:
	line = line.split("\t")
	try:
		genus2hightaxon_d[line[0]] = line[1]
	except IndexError:
		print(line, "not found")


#### Leaf2high taxon ####
#########################

for leaf in list(leavesfrompreds):
	tag = leaf.split("_")[0]
	if tag in taxarepl9:
		genus = taxarepl9[tag].split()[0]
	elif tag in genus2hightaxon_d:
		genus = tag
	else:
		print("WARNING: {} not present in genus-to-taxon translations. Leaves will be omitted from predictions.")
	try:
		curdict = {"hightaxon": genus2hightaxon_d[genus]} #this should not produce an exception
		preds_d[leaf].update(curdict)
	except KeyError:
		print("WARNING: {} not present in genus-to-taxon translations. Leaves will be omitted from predictions.")


#### Main ####
##############

#GROUP ASSIGNMENTS:
outgroups = {'Bacteria', 'Archaea'}
heterotrophs = {'Amoebozoa', 'Heterolobosea', 'Opisthokonta', 'SAR-Ciliophora', 'Parabasalia'}
opisthokonts = {'Metazoa', 'Fungi', 'Opisthokonta'} 	#Hectar / MultiLoc-animal
otherhetero = heterotrophs - opisthokonts 					#TargetP / MultiLoc-animal
primary = {'Rhodophyta', 'Viridiplantae', 'Glaucocystophyta'}	#TargetP / MultiLoc-plant
higherorder = {'SAR-Stramenopila', 'Haptophyta', 'Euglenozoa', 'SAR-Apicomplexa', 'SAR-Dinophyta', 'Cryptophyta', 'SAR-Rhizaria', }
eukaryote = heterotrophs | primary | higherorder
stramenopiles = {'SAR-Stramenopila'}		#Hectar / ASAFind
otherhigher = higherorder - stramenopiles 					#TargetP / ASAFind / MultiLoc-plant
ignored = {'unclassified', 'Undescribed'} | taxarepl9.keys()

missingset = set()
missingcount = 0

with open(prefix + '-preds_beta.txt','w') as outfilepredictions:
	for leaf in preds_d:
		fullgroup = preds_d[leaf]["hightaxon"]
		group = fullgroup.split("_")[0]
#######################################
#need this part
		if group in opisthokonts:
			prediction = ("{}: {}_//_{}".format(group, preds_d[leaf]["hectar"], preds_d[leaf]["ML2ANIMAL"]))
		elif group in otherhetero:
			prediction = ("{}: {}_//_{}".format(group, preds_d[leaf]["hectar"], preds_d[leaf]["ML2ANIMAL"]))
		elif group in primary:
			prediction = ("{}: {}_//_{}_//_{}".format(group, preds_d[leaf]["targetp"], preds_d[leaf]["signalp"], preds_d[leaf]["ML2PLANT"]))
		elif group in higherorder:
			prediction = ("{}: {}_//_{}_//_{}".format(group, preds_d[leaf]["hectar"], preds_d[leaf]["asafind"], preds_d[leaf]["ML2ANIMAL"]))
		elif group in otherhigher:
			prediction = ("{}: {}_//_{}_//_{}".format(group, preds_d[leaf]["hectar"], preds_d[leaf]["asafind"], preds_d[leaf]["ML2ANIMAL"]))
		else:
			#print(genus, "undetermined")
			prediction = "not determined"

########################################
		#prediction = re.sub('[():,]', '', prediction)
		#zde pridat predikci do taxa dictionary
		#outfileleaves.write("{}@{}__{}\n".format(query, fullgroup, prediction))
		outfilepredictions.write("{}\t{}\n".format(leaf, prediction))
	

print("prediction summaries written to file: {}-preds.txt".format(prefix))

#FINDING MISSING HIGHER TAXON DATA
if missingcount > 0:
	print("{} unassigned taxa found, consider rerunning script after updating the high_taxon_assignment.txt file with the following lines:".format(missingcount))
	with open("fetch_entrez_missing.txt", "w") as result:
		for item in missingset:
			result.write(item + "\n")

	from ete3 import NCBITaxa
	ncbi = NCBITaxa()
	
	missing = open("fetch_entrez_missing.txt").read().split("\n")
	name2taxid = ncbi.get_name_translator(missing)

	for genus in name2taxid:
		#retriev at least one species:
		descendants = ncbi.get_descendant_taxa(genus)
		lineage = ncbi.get_lineage(descendants[0])[2:7]
		names = ncbi.get_taxid_translator(lineage)
		rank = [names[taxid] for taxid in lineage]
		if "Eukaryota" in rank:
			rank.remove("Eukaryota")
		print("{}\t{}".format(genus, "_".join(rank)))

print("==============================================================")
print("Leaf renaming dictionary ready. Now to tree leaves renaming...")
print("Manual curation of {}-preds.txt needed. Please add 3rd column with possible preds:".format(prefix))
print("MT (mitochondrion), PT (plastid), CS (cytosol), dual (plastid+mitochondrion, amb (ambiguous), SP (signal peptide)")

if query_yes_no("Do you wish to continue?") is True:
	pass
else:
	quit("Algorithm interrupted manually.")

renaming = {}
with open("{}-preds.txt".format(prefix)) as f:
	for line in f:
		line = line.split("\t")
		print(line)
		if len(line) == 3:
			renaming[line[0]] = line[2]
		else:
			print("input improperly parsed")
#print(renaming)
quit("Algorithm not finished. Look up the ETE_branch_colours.py script")


#taxon renaming in the phylogenetic trees
inTrees = [s for s in os.listdir('.') if s.endswith('.tre') or s.endswith('.treefile')]
for currtree in inTrees:
	if currtree.endswith('-RXM.tre'):
		TREETYPE = "RAxML"
		if "bipartitions" in currtree:
			currtreename = currtree.split(".tre")[0].replace("RAxML_bipartitions.", "")
			#INVOKE TREE FUNCTION

			# tree_line = open(currtree).readline()
			# for key in taxa:
			# 	tree_line = tree_line.replace(key, taxa[key])
			# with open(currtreename + "-final.tre", "w") as result:
			# 	result.write(tree_line)
		else:
			print("skipping {}, not a tree file".format(currtree))

	elif currtree.endswith('.chain.con.tre'):
		TREETYPE = "PhyloBayes"
		currtreename = currtree.split(".tre")[0].replace(".chain.con", "")
		#INVOKE TREE FUNCTION

		# tree_line = open(currtree).readline()
		# for key in taxa:
		# 	tree_line = tree_line.replace(key, taxa[key])
		# with open(currtreename + "-PB-final.tre", "w") as result:
		# 	result.write(tree_line)

	elif currtree.endswith('.treefile'):
		TREETYPE = "IQtree"
		currtreename = currtree.split(".tre")[0].replace(".MaffTrimal.phy", "")
		#INVOKE TREE FUNCTION

		# tree_line = open(currtree).readline()
		# for key in taxa:
		# 	if key in tree_line:
		# 		tree_line = tree_line.replace(key, taxa[key])
		# 	elif key.replace("-","_") in tree_line:
		# 		modifkey = key.replace("-","_")
		# 		tree_line = tree_line.replace(modifkey, taxa[key])
		# 	elif key.replace("=","_") in tree_line:
		# 		modifkey = key.replace("=","_")
		# 		tree_line = tree_line.replace(modifkey, taxa[key])
		# 		#print(key, "found after modification")
		# with open(currtreename + "-IQT-final.tre", "w") as result:
		# 	result.write(tree_line)

	elif currtree.endswith('-ASA.tre'):
		TREETYPE = "AsaturA"
		currtreename = currtree.split(".tre")[0]
		#INVOKE TREE FUNCTION
		
		# tree_line = open(currtree).readline()
		# for key in taxa:
		# 	tree_line = tree_line.replace(key, taxa[key])
		# with open(currtreename + "-final.tre", "w") as result:
		# 	result.write(tree_line)

	else:
		print(currtree, " is a tree file of unknown origin (probably modified already?), skipping.")

print("Hotovo")
