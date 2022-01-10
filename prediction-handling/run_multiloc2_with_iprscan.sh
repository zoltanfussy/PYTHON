#run as: ~/localize/run_multiloc2_with_iprscan.sh some.fasta animal prefix
SCRATCHDIR="/tmp/zoli"
if [ $(grep -c "*" $1) -ne 0 ]; then
	echo "* character in input -- interpro will crash"
	echo "attempting to remove *, try again"
	sed -i 's/[*]\+//g' $1
	exit 1
else
	echo "input seems okay"
fi

interproscan.sh -T $SCRATCHDIR -appl Pfam -dp -i $1 -o interpro.out -f TSV -goterms -iprlookup
python /opt/multiloc2/git/MultiLoc2/MultiLoc2/src/multiloc2_prediction.py -fasta=$1 -origin=$2 -result=$3 -go=interpro.out -predictor=LowRes