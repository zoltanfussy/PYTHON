#./blastn.sh some.fasta
BLASTDIR="/opt/databases/nt" #2019-10-27
WDIR="/mnt/mokosz/home/zoli/DMND"

#to pass an infile argument
if [ $# -eq 1 ]
	then
		SAMPLE=$1
else
	SAMPLE="phylo.fasta"
fi

cd $WDIR

blastn -db $BLASTDIR/nt \
		-query $SAMPLE -out $SAMPLE.blast \
		-outfmt "6 qseqid staxids bitscore sseqid qcovs pident" \
		-max_target_seqs 5 -max_hsps 1 -evalue 1e-2 \
		-num_threads 4
