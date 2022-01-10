workdir='/mnt/mokosz/home/zoli/localize/'
outfile=${1%.*}'.ipsort.txt'
echo "output to $outfile"
cd $workdir

~/MitoFates/ipsort -type nonplant -i $1 -F -o $outfile
