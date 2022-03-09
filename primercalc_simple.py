def simple_tm(string):
	tm = 2*string.count("A") + 2*string.count("T") + 2*string.count("W") + \
		 4*string.count("G") + 4*string.count("C") + 4*string.count("S") - 5

	return tm


with open("primers.txt") as f:
	for l in f:
		data = l.strip().split("\t")
		primername = data[0]
		tm = simple_tm(data[1])
		print("primer {}\t Tm {} °C".format(primername, tm))
