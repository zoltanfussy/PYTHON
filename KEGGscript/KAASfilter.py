infile = open("query.ko").readlines()
out = open("filt.query.ko", "w")
n = 0
m = 0
for line in infile:
	n +=1
	if len(line.split()) == 2:
		print(line.split()[0])
		out.write(line)
		m += 1
print(n, m)
out.close()