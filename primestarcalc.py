with open("primers.txt") as f:
	for l in f:
		data = l.strip().split("\t")
		primername = data[0]
		seq = data[1]
		tm = 2*(seq.count("A") + seq.count("T")) + 4*(seq.count("C") + seq.count("G")) - 5
		print("primer {}\t Tm {} °C".format(primername, tm))
