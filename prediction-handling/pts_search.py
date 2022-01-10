import os,re,sys
from Bio import SeqIO

file = sys.argv[1]

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
elif os.path.isdir("/mnt/mokosz/home/zoli/"):
	home = "/mnt/mokosz/home/zoli/"
else:
	print("Please set a homedir")
defdir = ""
wd = home + defdir
#os.chdir(wd)

proteins = SeqIO.parse(file, 'fasta')

# #kinetoplastids
# pts1 = r'(S|A|G|C|N|P)(R|H|K|N|Q)(L|I|V|F|A|M|Y)'
# pts2 = r'^M\w{0,20}(R|K)(L|V|I)\w{5}(H|K|Q|R)(L|A|I|V|F|Y)'

#general
#pts1 = r'(S|A|C)(K|R|H|Q)(L|M)'
pts1 = r'[SAC][KRHQ][LM]'
#pts2 = r'^\w{1,21}R(L|I|V|Q)\w{2}(L|I|V|Q|H)(L|S|G|A)\w{1}(H|Q)(L|A)'
pts2 = r'^\w{1,21}R[LIVQ]\w{2}[LIVHQ][LSGA]\w{1}[HQ][LA]'

p1, p2, c = 0, 0, 0
with open(file.split(".")[0] + ".PTS.txt", 'w') as outlog, open(file.split(".")[0] + ".PTS.fasta", 'w') as out:
	for protein in proteins:
		c += 1
		if re.search(pts1, str(protein.seq)[-3:]):
			p1 += 1
			outlog.write("{}\tPTS1\n".format(protein.name))
			out.write('>{} @PTS1:{}\n{}\n'.format(protein.description, protein.seq[-3:], protein.seq))
		elif re.search(pts2, str(protein.seq)):
			p2 += 1
			match = re.search(pts2, str(protein.seq))
			outlog.write("{}\tPTS2\n".format(protein.name))
			out.write('>{} @PTS2:{}\n{}\n'.format(protein.description, match.group(), protein.seq))
		else:
			outlog.write("{}\tnone\n".format(protein.name))
		#	print(protein.description, '_____no PTS signal')

print("{} sequences processed".format(c))
print("Found {} PTS1 /{} PTS2 proteins".format(p1, p2))
