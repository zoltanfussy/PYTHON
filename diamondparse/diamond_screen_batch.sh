#/bin/bash
#PBS -N Diamond_nr_blastp
#PBS -q default
#PBS -l select=1:ncpus=8:mem=250gb:scratch_local=100gb
#PBS -l walltime=48:00:00
#PBS -m ae
#PBS -j oe
cat $PBS_NODEFILE

DATADIR="/storage/brno3-cerit/home/fussyz01/DMND"
#QUERY="DAPA_4.3.3.7.clust.fasta"
#add modules
module add diamond-0.8.29
module add python36-modules-gcc

#cp $DATADIR/refseq_ryby.dmnd $RESDIR/compared_annotations.fasta $SCRATCHDIR || exit 1
#cd $SCRATCHDIR
#diamond getseq --db nr.dmnd > nr.fasta
#./diamond_0.9.24 makedb --in nr.fasta -d nr_0.9.24
#cd $DATADIR
#cat refseq_ryby.fasta excavata.fasta > $SCRATCHDIR/refseq_ryby_excavata.fasta
#cp diamond_0.9.24 $SCRATCHDIR
#cd $SCRATCHDIR
#./diamond_0.9.24 makedb --in refseq_ryby_excavata.fasta -d refseq_excavata
#cp refseq_excavata.dmnd $DATADIR

cp $DATADIR/diamond_0.9.24 $SCRATCHDIR
cp $DATADIR/nr_0.9.24.dmnd $DATADIR/roots9.dmnd $DATADIR/MMETSPclust.dmnd $DATADIR/refseq_excavata.dmnd $SCRATCHDIR || exit 1
cp $DATADIR/*.fasta $SCRATCHDIR
cd $SCRATCHDIR

for QUERY in *fasta; do
	./diamond_0.9.24 blastp -p 8 -d roots9.dmnd -q $QUERY -o roots9.out -f 6 qseqid sseqid stitle evalue bitscore full_sseq --sensitive --max-target-seqs 15 --evalue 1e-5
	cp roots9.out $DATADIR
	./diamond_0.9.24 blastp -p 8 -d MMETSPclust.dmnd -q $QUERY -o mmetsp.out -f 6 qseqid sseqid stitle evalue bitscore full_sseq --sensitive --max-target-seqs 15 --evalue 1e-5
	cp mmetsp.out $DATADIR
	./diamond_0.9.24 blastp -p 8 -d refseq_excavata.dmnd -q $QUERY -o refseq.out -f 6 qseqid sseqid stitle evalue bitscore full_sseq --sensitive --max-target-seqs 15 --evalue 1e-5
	cp refseq.out $DATADIR
	./diamond_0.9.24 blastp -p 1 -d nr_0.9.24.dmnd -q $QUERY -o nr.out -f 6 qseqid sseqid stitle evalue bitscore full_sseq --sensitive --max-target-seqs 75 --evalue 1e-5
	cp nr.out $DATADIR

	cd $DATADIR
	python diamondparse.py --query $QUERY --results nr.out,refseq.out,mmetsp.out,roots9.out --pool True
	cd $SCRATCHDIR
done