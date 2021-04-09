from Bio import AlignIO

align = AlignIO.read("trans_multigene_v2_slowsites0.9.phy", "phylip-relaxed")
filt = align[:, :1] #has to be two ranges!
#print(filt)
seqs, length = len(align), align.get_alignment_length()
#print(align)
print(seqs,length)
for i,col in enumerate(range(1,length)):
	if i % 1000 == 0:
		print("processed", i)
	if align[:, col].count("X") < seqs:
		#print(align[:, col:col+1]) #just col won't work
		filt += align[:, col:col+1]

AlignIO.write(filt, "trans_multigene_v2_slowsites0.9_noX.phy", "phylip-relaxed")