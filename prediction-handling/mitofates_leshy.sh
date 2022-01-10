workdir='/mnt/mokosz/home/zoli/localize/'
outfile=${1%.*}'.mitofates.txt'
echo "output to $outfile"
cd $workdir
perl ~/MitoFates/MitoFates.pl $1 metazoa > $outfile
#choose from: fungi, metazoa or plant
