#./blastn.sh some.fasta
#BLASTDIR="/opt/databases/nt" #2019-10-27
BLASTDIR="/opt/databases/nt_auto/current/blast"
#WDIR="/mnt/mokosz/home/kika/rhizomastix_reassembly"
WDIR="/mnt/mokosz/home/zoli/DMND"

#to pass an infile argument
if [ $# -eq 1 ]
	then
		SAMPLE=$1
else
	SAMPLE="query.fasta"
fi

cd $WDIR

blastn -db $BLASTDIR/nt \
		-query $SAMPLE -out ${SAMPLE%.*}.blast \
		-outfmt "6 qseqid staxids bitscore sseqid qcovs pident" \
		-max_target_seqs 5 -max_hsps 1 -evalue 1e-10 \
		-num_threads 4
