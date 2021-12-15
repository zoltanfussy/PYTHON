
def simple_tm(string):
	tm = 2*string.count("A") + 2*string.count("T") + 2*string.count("W") + \
		 4*string.count("G") + 4*string.count("C") + 4*string.count("S") - 5

	return tm


primers = {">V9eug":
"ACCTTGTTACGACTTTTGC",
">V9FW2":
"TTTGTACACACCGCCC",
">V9meta":
"GTTACGACTTCTCSTTCCT",
">V9RV":
"CCTTCYGCAGGTTCACCTAC"}

for primer in primers:
	print("{}:\t{}".format(primer, simple_tm(primers[primer])))