from Bio import SeqIO,AlignIO
from numpy import log2

#FX

def Sobs(string, letter):
	lenstring = len(string)
	lettcount = string.count(letter)
	number = lettcount / lenstring
	if number != 0:
		value = -(number*log2(number))
		#value = number
	else:
		value = 0
	return value

def bitscorefreq(string, letter, weight):
	lenstring = len(string)
	lettcount = string.count(letter)
	number = lettcount / lenstring
	if number != 0:
		value = weight*number
	else:
		value = 0
	return value

#files:	
signalpfile = "plastid_refs-asafind.txt"
predisifile = "plastid_refs-predisi.txt"
predslfile = "plastid_refs-predsl.txt"

with open(signalpfile) as f:
	asafind = f.readlines()

with open(predisifile) as f:
	predisi = f.readlines()

with open(predslfile) as f:
	predsl = f.readlines()

#parse predictions
monstr_dic = {}
report_errors = False
print("Parsing signal peptide predictions and sequences into data dictionary...")
names = []
newnames = set()
for line in predisi:
	line = line.split("\t")
	#FASTA-ID	Cleavage Position	Signal Peptide ?	Score
	#[0]		[1]					[2]					[3]
	name = line[0]
	if line[0] not in names:
		names.append(name)
		if line[2] == "Y" and int(line[1]) < 101: #sometimes we see signal very far in the protein
			monstr_dic[name] = {"PrediSi site": int(line[1])+1, "PrediSi SPscore": float(line[3])}
		else:
			monstr_dic[name] = {"PrediSi site": "-", "PrediSi SPscore": float(line[3])}
	else:
		newnames.add(name)
		print("NEW NAME predisi! " + name)
print("->PrediSi: done")

names.sort()
#print(names)

for line in predsl:
	line = line.split("\t")
	#sequence id	mTP score	SP score	prediction	cleavage site
	#[0]			[1]			[2]			[3]			[4]			
	name = line[0]
	if name in names:
		#line[3] for nonplant prediction - complex algae
		if line[3] == "secreted":
			monstr_dic[name].update({"PredSL site": int(line[4])+1, "PredSL SPscore": float(line[2])})
		else:
			monstr_dic[name].update({"PredSL site": "-", "PredSL SPscore": float(line[2])})
	else:
		newnames.add(name)
		print("NEW NAME predsl! " + name)
print("->PredSL: done")

for line in asafind:
	#Identifier	SignalP	ASAfind cleavage position	ASAfind/SignalP cleavage site offset	ASAfind 20aa transit score	ASAfind Prediction
	#[0]		[1]		[2]							[3]										[4]							[5]	
	
	if not line.startswith('#') and len(line) != 0:
		line = line.split('\t')
		name = line[0]
		if name in names:
			try:
				sp = int(line[2])
				asa = sp + int(line[3])
			except ValueError:
				sp = "-"
				asa = "-"
			monstr_dic[name].update({'SignalP site': sp, 'ASAfind site': asa})
		else:
			newnames.add(name)
			print("NEW NAME asafind! " + name)
print("->ASAfind/SignalP: done")

if len(newnames) > 0:
	print("ERROR, new sequence names appeared!")
	report_errors = True

#load fastas:
fasta = SeqIO.parse("plastid_refs.fasta", "fasta")
seq_dic = {}
for seq in fasta:
	seq_dic[seq.name] = seq.seq
#print(seq_dic)


#process cleavage sites	
with open("plastid_signals.txt", "w") as out, open("cleavage_site.fasta", "w") as outfas:
	agreed = 0
	for name in names:
		#print(name)
		sequence = str(seq_dic[name])
		values = [monstr_dic[name]["PrediSi site"], monstr_dic[name]["PredSL site"], monstr_dic[name]["SignalP site"], monstr_dic[name]["ASAfind site"]]
		valueset = set(values)
		for i in valueset:
			if i != "-":
				if len(valueset - {"-"}) == 1 and values.count(i) == 2:
					agreed += 1
					print(name, i)
					outfas.write(">{}\n{}\n".format(name, sequence[int(i)-6:int(i)+19]))
				elif values.count(i) >=3:
					agreed += 1
					print(name, i)
					outfas.write(">{}\n{}\n".format(name, sequence[int(i)-6:int(i)+19]))
			
		out.write("{}\t{}\t{}\t{}\t{}\n".format(name, monstr_dic[name]["PrediSi site"], monstr_dic[name]["PredSL site"], 
			monstr_dic[name]["SignalP site"], monstr_dic[name]["ASAfind site"]))

	print(agreed, "sequences passed and written to alignment")

aacids = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X"]

alignmentcvel = AlignIO.read("cleavage_site_cvel.fasta", "fasta")
alignmentvbra = AlignIO.read("cleavage_site_vbra.fasta", "fasta")
#Rseq values calculated in Supplementary File S2 of the chromerid localizations paper
cvel_rseq = [0.628142594, 0.901165584, 1.679998503, 0.899323492, 2.865290913, 2.930248396, 1.218694579, 1.216631242, 0.723489751, 0.796375671, 0.62329697, 0.67311194, 0.660968002, 0.832345896, 0.930850003, 0.68832463, 0.701284019, 0.768918374, 0.590580854, 0.8231282, 0.710945491, 0.88885071, 0.651382747, 0.761732523, 0.743346288]
vbra_rseq = [0.817328158, 1.0450741, 1.916029505, 0.715175362, 3.370870006, 1.554297938, 0.863928534, 1.004913307, 0.97219529, 0.779193055, 0.870728277, 0.732232068, 0.787914214, 1.006489747, 0.892195734, 0.821817644, 0.999711101, 0.742618301, 0.895313802, 0.755506939, 0.848577712, 0.88178165, 0.706938911, 0.867174927, 0.857195538]

if cvel_rseq == []: #assuming Rseq values have not been calculated...
	with open("cvelmatrix.txt", "w") as cvelmatrix, open("vbramatrix.txt", "w") as vbramatrix:
		for aa in aacids:
			cveldata = [aa]
			vbradata = [aa]
			for i in range(25):
				#print(i)
				cveldata.append(Sobs(alignmentcvel[:, i], aa))
				vbradata.append(Sobs(alignmentvbra[:, i], aa))
			cvelstring = [str(x) for x in cveldata]
			vbrastring = [str(x) for x in vbradata]
			cvelmatrix.write("{}\n".format("\t".join(cvelstring)))
			vbramatrix.write("{}\n".format("\t".join(vbrastring)))
else:
	with open("cvelbitscorematrix.txt", "w") as cvelmatrix, open("vbrabitscorematrix.txt", "w") as vbramatrix:
		for aa in aacids:
			cveldata = [aa]
			vbradata = [aa]
			for i in range(25):
				#print(i)
				cveldata.append(bitscorefreq(alignmentcvel[:, i], aa, cvel_rseq[i]))
				vbradata.append(bitscorefreq(alignmentvbra[:, i], aa, vbra_rseq[i]))
			cvelstring = [str(x) for x in cveldata]
			vbrastring = [str(x) for x in vbradata]
			cvelmatrix.write("{}\n".format("\t".join(cvelstring)))
			vbramatrix.write("{}\n".format("\t".join(vbrastring)))


