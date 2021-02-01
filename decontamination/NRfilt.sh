#./blastn.sh some.fasta
BLASTDIR="/opt/databases/nr_auto/current/diamond" #2019-10-27
WDIR="/mnt/mokosz/home/zoli/DMND"

#to pass an infile argument
if [ $# -eq 1 ]
	then
		SAMPLE=$1
else
	SAMPLE="phylo.fasta"
fi

cd $WDIR

#diamond blastx
diamond blastx -d $BLASTDIR/nr.dmnd \
	-q $SAMPLE -o $SAMPLE.dmnd.out \
	-f 6 qseqid bitscore sseqid qcovhsp pident qlen length --sensitive \
	--max-target-seqs 1 --evalue 1e-5 \
	-p 10

diamond blastx -d roots9.dmnd \
	-q $SAMPLE -o $SAMPLE.roots9.dmnd.out \
	-f 6 qseqid bitscore sseqid qcovhsp pident qlen length --sensitive \
	--max-target-seqs 1 --evalue 1e-5 \
	-p 10

#now to blobtools
python taxify_DMND_nr_gz.py -i $SAMPLE.dmnd.out
python taxify_DMND_nr_gz.py -i $SAMPLE.roots9.dmnd.out
#blast output already in taxified format:
