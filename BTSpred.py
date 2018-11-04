print("Suggested usage: python BTSpred.py -in INFILE")
print("make sure predisi.txt and predsl.txt files are ready")


import argparse
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#os.system('source ~/.bash_profile')

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-in', '--infile', help='fasta to be processed', default="chromera.fasta")

args = parser.parse_args()

filename = args.infile.split(".")[0]

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
		raise ValueError("invalid default answer: '%s'" % default)

	while True:
		sys.stdout.write(question + prompt)
		if sys.platform in ["darwin", "win32"]:
			choice = input().lower()
		elif sys.platform.startswith("linux"):
			choice = raw_input().lower()
		else:
			quit("unrecognized OS, check how to use raw input")
		
		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "
							 "(or 'y' or 'n').\n")


monstr_dic = {}
###########################################################################
#		read in seq names to extract nt/aa models of interest:		   #
###########################################################################
inproteins = SeqIO.parse(args.infile, 'fasta')
names = []
for sequence in inproteins:
	annotation = " ".join(sequence.description.split(" ")[1:])
	name = sequence.name
	if "-reversed" in name:
		name = name.replace("-reversed", "")
	if name not in names:	
		names.append(name)
	seq = sequence.seq
	monstr_dic[name] = 	{"Annotation": annotation, "Name": name, "Seq": seq}


#print("Now printing control table, number of items per dictionary key...")
#for sequence in monstr_dic.keys():
#	print(sequence, len(monstr_dic[sequence]))

###########################################################################
#		  read in predictions from PrediSi, PredSL, TMHMM:			   #
###########################################################################
predisi = open("predisi.txt").readlines()
predsl = open("predsl.txt").readlines()


print("parsing PrediSi and PredSL predictions and sequences into data dictionary")
for line in predisi:
	line = line.split("\t")
	#FASTA-ID	Cleavage Position	Signal Peptide ?	Score
	#[0]		[1]					[2]					[3]
	if line[0] in names:
		name = line[0]
		if line[2] == "Y" and int(line[1]) < 101: #sometimes we see signal very far in the protein
			monstr_dic[name].update({"PrediSi site": line[1], "PrediSi SPscore": float(line[3])})
		else:
			monstr_dic[name].update({"PrediSi site": "-", "PrediSi SPscore": float(line[3])})

for line in predsl:
	line = line.split()
	#sequence id	mTP score	SP score	prediction	cleavage site
	#[0]			[1]			[2]			[3]			[4]			
	if line[0] in names:
		name = line[0]
		#line[3] for nonplant prediction - complex algae
		if line[3] == "secreted":
			monstr_dic[name].update({"PredSL site": line[4], "PredSL SPscore": float(line[2])})
		else:
			monstr_dic[name].update({"PredSL site": "-", "PredSL SPscore": float(line[2])})


###########################################################################
#					  transit peptide detection:						 #
###########################################################################

print("preparing output for MultiLoc2")

with open('cleaved-' + args.infile, 'a') as result:
	for protein in monstr_dic.keys():
		#print(protein, monstr_dic[protein])
		try:
			print(protein)
			SP_si = monstr_dic[protein]["PrediSi site"]
			SP_sl = monstr_dic[protein]["PredSL site"]
			if SP_si == "-":
				if SP_sl != "-":
					SP_sl = int(SP_sl)
					SP_cleaved_sl = '>{}\n{}\n'.format(protein, monstr_dic[protein]["Seq"][SP_sl:])
					#SP_cleaved_sl = '>{}@PredSL\n{}\n'.format(protein, monstr_dic[protein]["Seq"][SP_sl:])
					result.write(SP_cleaved_sl)
			elif SP_sl == "-":
				SP_si = int(SP_si)
				SP_cleaved_si = '>{}\n{}\n'.format(protein, monstr_dic[protein]["Seq"][SP_si:])
				#SP_cleaved_si = '>{}@PrediSi\n{}\n'.format(protein, monstr_dic[protein]["Seq"][SP_si:])
				result.write(SP_cleaved_si)
			elif int(SP_si) == int(SP_sl):
				SP_si = int(SP_si)
				SP_sl = int(SP_sl)
				SP_cleaved = '>{}\n{}\n'.format(protein, monstr_dic[protein]["Seq"][SP_si:])
				result.write(SP_cleaved)
			else:
				SP_si = int(SP_si)
				SP_sl = int(SP_sl)
				SP_cleaved_si = '>{}\n{}\n'.format(protein, monstr_dic[protein]["Seq"][SP_si:])
				SP_cleaved_sl = '>{}@PredSL\n{}\n'.format(protein, monstr_dic[protein]["Seq"][SP_sl:])
				result.write(SP_cleaved_si)
				result.write(SP_cleaved_sl)
		except ValueError as VE:
			with open('result_errors.txt', 'a') as errors:
				errors.write('SP cleavage error: {}\n{}\n\n'.format(monstr_dic[protein][5], str(VE)))

print("cleaved proteins written to {}... now running MultiLoc2 plant option to determine transit peptide probabilities.".format("cleaved-" + args.infile))

if query_yes_no("Do you wish to perform MultiLoc2 on cleaved sequences?") is True:
	os.system('python /home/manager/MultiLoc2/src/multiloc2_prediction.py -fasta={} -predictor=LowRes -origin=plant -result=pred-ML2-cleave.txt -output=simple'.format("cleaved-" + args.infile))
#use the following if your script failed, but you have made the new MultiLoc2 predictions
#if os.path.isfile("pred-ML2-cleave.txt") == False:
#	os.system('python /home/manager/MultiLoc2/src/multiloc2_prediction.py -fasta={} -predictor=LowRes -origin=plant -result=pred-ML2-cleave.txt -output=simple'.format("cleaved-" + args.infile))

print("MultiLoc2 finished, now adding predictions to data dictionary")

with open("pred-ML2-cleave.txt") as infileML2:
	ML2C = infileML2.read().split('\n')

alternativecleavage = open("alternativecleavage-" + filename + ".txt", "w")
alternativecleavage_dict = {}

for line in ML2C:
	if len(line.split('\t')) != 1:
		line = line.split('\t')
		preds = sorted(line[1:])
		name = line[0]
		#['chloroplast: ', 'cytoplasmic: ', 'mitochondrial: ', 'nuclear: ', 'secretory pathway: ']
		if "@PredSL" in name or "@PrediSi" in name:
			alternativecleavage.write("alternative cleavage for: {}\n".format("\t".join(line)))
			name = name.split("@")[0]
			try:
				cTP = preds[0].split()[1]
				mTP = preds[2].split()[1]
				combinedMLpred = float(cTP) + float(mTP)
				if combinedMLpred > 0.3:
					ourguess = "BIPARTITE, plastid"
				else:
					ourguess = "secretory"
				#print(cTP, mTP)
				if name in alternativecleavage_dict: 
					if combinedMLpred > alternativecleavage_dict[name]: #here we choose the "better prediction"
						monstr_dic[name].update({"ML2post-cleavage": "pt: {} (mt: {})".format(cTP, mTP), "ourguess": ourguess, "ourscore": combinedMLpred})
				else:
					alternativecleavage_dict[name] = combinedMLpred
					monstr_dic[name].update({"ML2post-cleavage": "pt: {} (mt: {})".format(cTP, mTP), "ourguess": ourguess, "ourscore": combinedMLpred})

			except ValueError as VE:
				with open('result_errors.txt', 'a') as errors:
					errors.write('ML2 processing error: {}\n{}\n\n'.format(line, str(VE)))

		else:
			try:
				cTP = preds[0].split()[1]
				mTP = preds[2].split()[1]
				combinedMLpred = float(cTP) + float(mTP)
				if combinedMLpred > 0.3:
					ourguess = "BIPARTITE, plastid"
				else:
					ourguess = "secretory"
				#print(cTP, mTP)
				monstr_dic[name].update({"ML2post-cleavage": "pt: {} (mt: {})".format(cTP, mTP), "ourguess": ourguess, "ourscore": combinedMLpred})
			except ValueError as VE:
				with open('result_errors.txt', 'a') as errors:
					errors.write('ML2 processing error: {}\n{}\n\n'.format(line, str(VE)))
print("Transit peptide predictions made, check alternativecleavage" + filename + ".txt for ambiguous predictions.")
alternativecleavage.close()


print("Now writing results...")
outtable = 'FINALTABLE-' + filename + '.txt'

#solution without pandas
with open('FINALTABLE-' + filename, 'w') as outtable, open('FINALPREDS-' + filename, 'w') as outpreds:
	for name in monstr_dic:
		item = monstr_dic[name]
		ML2pred = item.get("ML2post-cleavage", "-")
		ourplastidpred = item.get("ourguess", "non-plastid")
		ourplastidscore = "%.3f" % ((max(item["PrediSi SPscore"], item["PredSL SPscore"]) + 3*item.get("ourscore", 0))/4)
		outtable.write("{}\t{}\t{}\t{}\n".format(item["Name"],item["PrediSi site"],item["PredSL site"],ML2pred))
		outpreds.write("{}\t{}\t{}\n".format(item["Name"],ourplastidpred, ourplastidscore))
		print(item["Name"], item["PrediSi SPscore"], item["PredSL SPscore"], item.get("ourscore", 0))


print("Writing done. Hooray.")
