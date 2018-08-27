from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt

scaffolds = SeqIO.parse("./seq-extend-BLASTXML/CryptoDB-36_CveliaCCMP2878_Genome.fasta", "fasta")

with open("candidate-centromeres.fasta", "w") as result:
	for rec in scaffolds:
		sequence = rec.seq
		if len(sequence) > 200000:
			plotcreate = True
			plotdata = {}
		else:
			plotcreate = False
		threek_windows = [sequence[start:start+3000] for start in range(0, len(sequence), 3000)]
		windowcount = 0
		for threekindex,i in enumerate(threek_windows):
			lowgc_windows_num = 0
			hundred_windows = [i[start:start+100] for start in range(0, len(i), 100)]
			for hunindex,window in enumerate(hundred_windows):
				window = str(window).replace("N","")
				GCperc = GC(window)
				if plotcreate:
					windowID = (threekindex)*30 + hunindex + 1
					plotdata[windowID] = GCperc
				if len(window) > 20 and GCperc < 32: #can be 32% as in the Diner et al 2017 paper
					lowgc_windows_num += 1
					#print("{} {}".format(windowID,window))
			if lowgc_windows_num > 10:
				result.write(">{}_candidate centromere {}-{}\n{}\n".format(rec.name, (threekindex)*3000, (threekindex+1)*3000, rec.seq))
		# if plotcreate: #only runs in jupyter - graphical output
			# plotvalues = [int(plotdata[i]) for i in range(1, len(plotdata.keys()))]
			# plt.plot(range(1, len(plotdata.keys())), plotvalues)
			# plt.show()
