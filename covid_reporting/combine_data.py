import tarfile
import io,os
from Bio import SeqIO

#5000 sequences 
#		 	-2021.07.07
# 2021.08.07-2021.10.14
# 2021.10.15-2021.11.21
# 2021.11.22-2021.12.23
# 2021.12.24-2022.01.19
# 2022.01.20-2022.02.15
# 2022.02.16-2022.03.

#set local computer here!
local = ["zoliq", "morpholino", "faculty"][0]
if local == "zoliq":
	os.chdir("/Volumes/zoliq data/Downloads/covid_seq")
else:
	os.chdir("~/Downloads/covid_seq")

with open("gisaid_combined.fasta", "wt") as result:
	#just reset this file at each run, below I am appending from multiple sources
	pass

metadata_dict = {}
headers_dict = {}
dupes = set()
gisaiddupes = set()
writtendupes = set()
known_dupes = set(["hCoV-19/CzechRepublic/UMTM270158/2021", 
				   "hCoV-19/CzechRepublic/UMTM228154/2021", 
				   "hCoV-19/CzechRepublic/NRL_14022/2021",
				   "hCoV-19/CzechRepublic/NRL_13996/2021",
				   "hCoV-19/CzechRepublic/NRL_14016/2021",
				   "hCoV-19/CzechRepublic/NRL_14025/2021",
				   "hCoV-19/CzechRepublic/NRL_13991/2021",
				   "hCoV-19/CzechRepublic/NRL_14029/2021",
				   "hCoV-19/CzechRepublic/NRL_13989/2021",
				   "hCoV-19/CzechRepublic/NRL_13995/2021",
				   "hCoV-19/CzechRepublic/NRL_14010/2021",
				   "hCoV-19/CzechRepublic/13SVAN94/2021",
				   "hCoV-19/CzechRepublic/UMTM335807/2021",
				   "hCoV-19/CzechRepublic/UMTM331673/2021",
				   "hCoV-19/CzechRepublic/UMTM331796/2021",
				   "hCoV-19/CzechRepublic/UMTM335807/2021",
				   "hCoV-19/CzechRepublic/NRL_12112/2021",
				   ])
#fasta_dict = {} >> this is probably unnecessary and I will append sequences to a result file

"""
#this part for extracted data
gisaid_dirs = [x for x in os.listdir(".") if x.startswith("gisaid_auspice") and not x.endswith("tar")]
print("Merging from:\n{}".format(gisaid_dirs))

for d in gisaid_dirs:
	metadata = [x for x in os.listdir(d) if x.endswith("tsv")][0]
	fasta = [x for x in os.listdir(d) if x.endswith("fasta")][0]
	with open("{}/{}".format(d, metadata)) as f:
		for l in f:
			seqid = l.split("\t")[0]
			if len(seqid) == 0:
				continue
			metadata_dict[seqid] = l
			l = l.split("\t")
			headers_dict[seqid] = ">{}|{}|{}".format(seqid, l[2], l[4])
	with open("gisaid_combined.fasta", "at") as result:
		for seq in SeqIO.parse("{}/{}".format(d, fasta), "fasta"):
			if seq.name in headers_dict:
				result.write("{}\n{}\n".format(headers_dict[seq.name], seq.seq))

#finally print metadata too:
with open("gisaid_combined.tsv", "wt") as result:
	for key in metadata_dict:
		result.write(metadata_dict[key])
"""

#this part for tar archives:
gisaid_dirs = [x for x in os.listdir(".") if x.startswith("gisaid_auspice") and x.endswith("tar")]
print("Merging from:\n{}".format(gisaid_dirs))

for d in gisaid_dirs:
	print(d)
	with tarfile.open(d) as tar:
		metadata = [x for x in tar.getmembers() if x.name.endswith("tsv")]
		fasta = [x for x in tar.getmembers() if x.name.endswith("fasta")]

		with tar.extractfile(metadata[0]) as f:
			content = io.TextIOWrapper(f)
			for l in content:
				seqid = l.split("\t")[0]
				if len(seqid) == 0:
					continue
				if seqid.startswith("hCoV-19") or seqid == "strain":
					data = l.split("\t")
					gisaid = data[2]
					seqidlen = "{}|{}".format(l.split("\t")[0], l.split("\t")[13])
					if seqid in headers_dict and seqid != "strain":
						print("Warning, duplicate seqID", seqid)
						dupes.add(seqid)
						gisaiddupes.add(gisaid)
					metadata_dict[gisaid] = l
					headers_dict[seqid] = ">{}|{}|{}".format(seqid, gisaid, data[4])
					headers_dict[seqidlen] = ">{}|{}|{}".format(seqid, gisaid, data[4])
					last_seqid = data[2]
				else:
					#print(last_seqid, l)
					metadata_dict[last_seqid] = "{}{}".format(metadata_dict[last_seqid].replace("\n", " "), l)
		with open("gisaid_combined.fasta", "at") as result,\
			tar.extractfile(fasta[0]) as f:
			content = io.TextIOWrapper(f)
			for seq in SeqIO.parse(content, "fasta"):
				if seq.name in dupes:
					seqidlen = "{}|{}".format(seq.name, len(seq.seq))
					if seqidlen not in writtendupes:
						result.write("{}\n{}\n".format(headers_dict[seqidlen], seq.seq))
						writtendupes.add(seqidlen)
						writtendupes.add(headers_dict[seqidlen].split("|")[1])
				elif seq.name in headers_dict:
					result.write("{}\n{}\n".format(headers_dict[seq.name], seq.seq))
				else:
					print("{} not in headers_dict".format(seq.name))

print("Known duplicates: {}".format(",".join([x for x in known_dupes])))

#finally print metadata too:
writtendupes = set()

with open("gisaid_combined.tsv", "wt") as result:
	for key in metadata_dict:
		if key in gisaiddupes:
			if key not in writtendupes:
				result.write(metadata_dict[key])
				writtendupes.add(key)
		else:
			result.write(metadata_dict[key])


#now xz the result
print("Please wait for compression to finish...")
if local == "zoliq":
	os.system("/Users/zoliq/anaconda3/bin/xz -f gisaid_combined.fasta") #on zoliq
elif local == "morpholino":
	os.system("/opt/anaconda3/bin/xz -f gisaid_combined.fasta") #on morpholino
elif local == "faculty":
	os.system("/Users/morpholino/opt/anaconda3/bin/xz -f gisaid_combined.fasta") #on faculty mac
print("Done!")
