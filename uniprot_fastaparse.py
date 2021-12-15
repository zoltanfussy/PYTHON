import re

orgnpattern = r"OS=(.*)OX="

with open("some.fasta") as f:
	result = f.read()
	#print(result)
rawresult = result.replace("tr|", "").replace("sp|", "")
resultstring = ""
for line in rawresult.split(">")[1:]:
	try:
		orgn = re.search(orgnpattern, line).group(1)
		orgn = "_".join(orgn.split()[:2])
		acc_entry = line.split()[0].replace("|", "@")
		annot = " ".join(line.split(" ")[1:]).split("OS=")[0]
		seq = "".join(line.split("\n")[1:])
		resultstring += f">{orgn}_{acc_entry} {annot}\n{seq}\n"
	except IndexError:
		pass

with open("some_converted.fasta", "w") as out:
	out.write(resultstring)