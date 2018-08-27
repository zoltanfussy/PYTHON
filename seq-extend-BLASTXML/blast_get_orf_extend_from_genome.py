#!/usr/bin/python3
#!!! Check parsing record.query in the blast_parser function (5x) !!!
from Bio import SeqIO
import re
from Bio.Blast import NCBIXML

#owncloudpath = "/Users/zoliq/ownCloud"
owncloudpath = "/Volumes/zoliq data/ownCloud"
fasta = SeqIO.parse(owncloudpath + '/large/seq-extend-BLASTXML/CryptoDB-36_VbrassicaformisCCMP3155_Genome.fasta', 'fasta')
queryfasta = SeqIO.parse(owncloudpath + '/large/seq-extend-BLASTXML/vitrella_refs.fasta', 'fasta')
nt_out = open(owncloudpath + '/large/seq-extend-BLASTXML/vitrella_refs-elong_nt_genome.txt', 'w')
aa_out = open(owncloudpath + '/large/seq-extend-BLASTXML/vitrella_refs-elong_aa_genome.txt', 'w')
err_out = open(owncloudpath + '/large/seq-extend-BLASTXML/vitrella_refs_errors_genome.txt', 'w')
logfile = open(owncloudpath + '/large/seq-extend-BLASTXML/vitrella_refs_logfile_genome.txt', 'w')
result_handle = open(owncloudpath + '/large/seq-extend-BLASTXML/vbra_refs_blast_genome.xml')
blast_records = NCBIXML.parse(result_handle)

maxmismatches = 11

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'B', 'TAG':'B',
    'TGC':'C', 'TGT':'C', 'TGA':'B', 'TGG':'W'}

def translation(sequence):
    cut_seq = []
    for i in range(0,len(sequence)-2,3):
        cut_seq.append(sequence[i:i+3])
    aa = []
    for codon in cut_seq:
        if 'N' in codon:
            aa.append('X')
        else:
            aa.append(gencode[codon])
    return ''.join(aa)

def blast_parser(blast_records):
	
	result = {}
	errors = []
	for record in blast_records:
		hspcount = 0
		try:
			logfile.write(record.query + "\n")
			best = record.alignments[0]
			hit2query[best.hit_id] = record.query
			min_sstart = False
			max_send = False
			min_qstart = False
			max_qend = False
			frame = best.hsps[0].frame[1]
			if best.hsps[0].expect > 0.01:
				err_out.write('{}:\ttoo high evalue\n'.format(record.query.split('Monsieur Guillotin')[0]))
			else:
				for hsp in best.hsps:
					hspcount += 1
					logfile.write("hsp{}\n".format(hspcount))
					mismatches = hsp.align_length - (hsp.gaps  + hsp.identities)
					if mismatches < maxmismatches:
						if frame == hsp.frame[1]:
							if not min_qstart:
								min_qstart = hsp.query_start
								if frame in [1, 2, 3]:
									min_sstart = best.hsps[0].sbjct_start
								else:
									min_sstart = best.hsps[0].sbjct_end
							if not max_qend:
								max_qend = hsp.query_end
								if frame in [1, 2, 3]:
									max_send = best.hsps[0].sbjct_end
								else:
									max_send = best.hsps[0].sbjct_start
							if min_qstart > hsp.query_start:
								min_qstart = hsp.query_start
								if frame in [1, 2, 3]:
									min_sstart = hsp.sbjct_start
								else:
									min_sstart = hsp.sbjct_end
							if max_qend < hsp.query_end:
								max_qend = hsp.query_end
								if frame in [1, 2, 3]:
									max_send = hsp.sbjct_end
								else:
									max_send = hsp.sbjct_start
						else:
							errors.append(record.query.split('Monsieur Guillotin')[0])
							#find out if this is the leftmost hsp!
							
							if hsp.frame[1] in [1, 2, 3]:
								if min_sstart > hsp.sbjct_start:
									logfile.write("FW strand, good HSP upstream from the last good HSP ({} srv. {})\n".format(hsp.sbjct_start, min_sstart))
									logfile.write("changing frame from {} to {}\n".format(frame, hsp.frame[1]))
									frame = hsp.frame[1]	
							else:
								if min_sstart < hsp.sbjct_end:
									logfile.write("RV strand, good HSP upstream from the last good HSP ({} srv. {})\n".format(hsp.sbjct_end, min_sstart))
									logfile.write("changing frame from {} to {}\n".format(frame, hsp.frame[1]))
									frame = hsp.frame[1]

							if min_qstart > hsp.query_start:
								min_qstart = hsp.query_start
								if frame in [1, 2, 3]:
									min_sstart = hsp.sbjct_start
								else:
									min_sstart = hsp.sbjct_end
							if max_qend < hsp.query_end:
								max_qend = hsp.query_end
								if frame in [1, 2, 3]:
									max_send = hsp.sbjct_end
								else:
									max_send = hsp.sbjct_start
					else:
						logfile.write("too many mismatches: {}\n".format(mismatches))
						#print("{}\n{}\n{}\n".format(hsp.query, hsp.match, hsp.sbjct)) #you may print the alignments

				if frame in [1, 2, 3]:
					result[record.query.split('Monsieur Guillotin')[0]] = [min_sstart, max_send, frame, best.hit_id, 
						record.query_length, min_qstart, max_qend]
				else:
					result[record.query.split('Monsieur Guillotin')[0]] = [max_send, min_sstart, frame, best.hit_id, 
						record.query_length, min_qstart, max_qend]
		except: #stupid for debugging, I suggest using KeyError
			err_out.write('{}:\tno hit found\n'.format(record.query.split('Monsieur Guillotin')[0]))
	errors = set(errors)
	for i in errors:
		err_out.write('{}:\thsps frames do not correspond\n'.format(i))
	return result

hit2query = {}
blast_dict = blast_parser(blast_records)
#print(blast_dict)
queryseq_dict = {}
for s in queryfasta:
	queryseq_dict[s.name] = (s.seq, len(s.seq))

for contig in fasta:
	if contig.name in hit2query:
		key = contig.name
		queryid = hit2query[contig.name]
		value = blast_dict[queryid]
		frame = value[2]
		ref_name = queryid
		min_qstart_nt = value[5]*3
		if frame in [1, 2, 3]:
			print(contig.name + '_____forward')
			seq_start = value[0]-1
			seq_end = value[1]
			prev_stop = seq_start - 3
			while translation(contig.seq[prev_stop:prev_stop+3]) != 'B': #B denotes stopcodon, see gencode table at the beginning
				if prev_stop > 2:
					prev_stop = prev_stop - 3
				else:
					prev_stop = prev_stop
					break
			else:
				prev_stop = prev_stop
			if 'M' not in translation(contig.seq[prev_stop:seq_start-1]):
				new_start = seq_start #we dont want incomplete N termini
				# if translation(contig.seq[seq_start:seq_start+3]) == 'M':
				# 	new_start = seq_start
				# else:
				# 	new_start = prev_stop + 3
			else:
				new_start = prev_stop + 3*translation(contig.seq[prev_stop:seq_start-1]).find('M')
			if translation(contig.seq[seq_end:seq_end+3]) == 'B':
				seq_end = seq_end
			else:
				while translation(contig.seq[seq_end:seq_end+3]) != 'B':
					if seq_end < len(contig.seq) - 3:
						seq_end = seq_end + 3
					else:
						seq_end = seq_end
						break
				else:
					seq_end = seq_end
			nucleotides = contig.seq[new_start:seq_end+3]
			N_extension = contig.seq[new_start:max(seq_start-min_qstart_nt+3, new_start)] #we want some overhang, but not reverse sequence, +3 necessary to keep all aa
			#protein = translation(nucleotides)[:-1] #original extraction, good for sequences without introns
			protein = translation(N_extension) + "^" + queryseq_dict[queryid][0]
			if len(protein) - 1 > queryseq_dict[queryid][1]:
				print("found longer sequence, writing")
				nt_out.write('>{}__{}\n{}\n'.format(contig.name, ref_name, nucleotides))
				aa_out.write('>{}__{}\n{}\n'.format(contig.name, ref_name, protein))
			else:
				print("hit shorter than query")
		else:
			print(contig.name + '_____reverse')
			reverse = contig.seq.reverse_complement()
			seq_start = len(reverse) - value[1]
			seq_end = len(reverse) - value[0] + 1
			prev_stop = seq_start - 3
			
			while translation(reverse[prev_stop:prev_stop+3]) != 'B':
				if prev_stop > 2:
					prev_stop = prev_stop - 3
				else:
					prev_stop = prev_stop
					break
			else:
				prev_stop = prev_stop
			if 'M' not in translation(reverse[prev_stop:seq_start-1]):
				new_start = seq_start #we dont want incomplete N termini
				# if translation(reverse[seq_start:seq_start+3]) == 'M':
				# 	new_start = seq_start
				# else:
				# 	new_start = prev_stop + 3
			else:
				new_start = prev_stop + 3*translation(reverse[prev_stop:seq_start-1]).find('M')
			if translation(reverse[seq_end:seq_end+3]) == 'B':
				seq_end = seq_end
			else:
				while translation(reverse[seq_end:seq_end+3]) != 'B':
					if seq_end < len(reverse) - 3:
						seq_end = seq_end + 3
					else:
						seq_end = seq_end
						break
				else:
					seq_end = seq_end
			nucleotides = reverse[new_start:seq_end+3]
			N_extension = reverse[new_start:max(seq_start-min_qstart_nt+3, new_start)] #we want some overhang, but not reverse sequence, +3 necessary to keep all aa
			#protein = translation(nucleotides)[:-1] #original extraction, good for sequences without introns
			protein = translation(N_extension) + "^" + queryseq_dict[queryid][0]
			if len(protein) - 1 > queryseq_dict[queryid][1]:
				print("found longer sequence, writing")
				nt_out.write('>{}__{}\n{}\n'.format(contig.name, ref_name, nucleotides))
				aa_out.write('>{}__{}\n{}\n'.format(contig.name, ref_name, protein))
			else:
				print("hit shorter than query")
	else:
		pass

nt_out.close()
aa_out.close()