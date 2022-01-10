targetp -fasta $1 -org non-pl -format short #-mature
outfile=${1%.*}
#~/localize/targetp.sh file
#-mature <FILE>
#-batch int
mv "$outfile"_summary.targetp2 "$outfile".targetp2.txt