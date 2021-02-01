#/bin/bash
#PBS -N Diamond_nr_blastp
#PBS -q default
#PBS -l select=1:ncpus=8:mem=250gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -j oe
cat $PBS_NODEFILE

DATADIR="/storage/brno3-cerit/home/fussyz01/DMND"
#qsub -v SAMPLE=hspN.fasta diamond_screen_one.sh
#SAMPLE="whatever.fasta"
#add modules
module add diamond-0.8.29
module add python36-modules-gcc

#cp $DATADIR/refseq_ryby.dmnd $RESDIR/compared_annotations.fasta $SCRATCHDIR || exit 1
#cd $SCRATCHDIR
#diamond getseq --db nr.dmnd > nr.fasta
#./diamond_0.9.24 makedb --in nr.fasta -d nr_0.9.24
cd $DATADIR
cat refseq_ryby.fasta hapto_clust.fasta > $SCRATCHDIR/refseq_hapto.fasta
cp diamond_0.9.24 $SCRATCHDIR
cd $SCRATCHDIR
./diamond_0.9.24 makedb --in refseq_hapto.fasta -d refseq_hapto
cp refseq_hapto.dmnd $DATADIR

cp $DATADIR/diamond_0.9.24 $SCRATCHDIR
#cp $DATADIR/refseq_hapto.dmnd $SCRATCHDIR || exit 1
cp $DATADIR/nr_0.9.24.dmnd $DATADIR/roots9.dmnd $DATADIR/MMETSPclust.dmnd $SCRATCHDIR || exit 1
cp $DATADIR/phaant_ref_aa.dmnd $DATADIR/phaglo_ref_aa.dmnd $SCRATCHDIR || exit 1
cp $DATADIR/$SAMPLE $SCRATCHDIR
cd $SCRATCHDIR

#smaller to bigged databases...
./diamond_0.9.24 blastp -p 10 -d roots9.dmnd -q $SAMPLE -o roots9.out -f 6 qseqid sseqid evalue ppos stitle full_sseq --sensitive --max-target-seqs 80 --evalue 1e-5
cp roots9.out $DATADIR
./diamond_0.9.24 blastp -p 10 -d refseq_hapto.dmnd -q $SAMPLE -o refseq.out -f 6 qseqid sseqid evalue ppos stitle full_sseq --sensitive --max-target-seqs 80 --evalue 1e-5
cp refseq.out $DATADIR
./diamond_0.9.24 blastp -p 10 -d phaant_ref_aa.dmnd -q $SAMPLE -o phaant.out -f 6 qseqid sseqid evalue ppos stitle full_sseq --sensitive --max-target-seqs 15 --evalue 1e-5
cp phaant.out $DATADIR
./diamond_0.9.24 blastp -p 10 -d phaglo_ref_aa.dmnd -q $SAMPLE -o phaglo.out -f 6 qseqid sseqid evalue ppos stitle full_sseq --sensitive --max-target-seqs 15 --evalue 1e-5
cp phaglo.out $DATADIR
#./diamond_0.9.24 blastp -p 10 -d MMETSPclust.dmnd -q $SAMPLE -o mmetsp.out -f 6 qseqid sseqid evalue ppos stitle full_sseq --sensitive --max-target-seqs 30 --evalue 1e-5
#cp mmetsp.out $DATADIR
./diamond_0.9.24 blastp -p 10 -d nr_0.9.24.dmnd -q $SAMPLE -o nr.out -f 6 qseqid sseqid evalue ppos stitle full_sseq --sensitive --max-target-seqs 80 --evalue 1e-5
cp nr.out $DATADIR

cd $DATADIR
python3 diamondparse.py --query $SAMPLE --pool True --skip_dupes --results nr.out,refseq.out,roots9.out,phaant.out,phaglo.out || break
#python diamondparse_server.py --query $SAMPLE --pool True --results phylodb.out || break
#originally --results genescreen - parses all the above outputs
