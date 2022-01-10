#use _allpreds.sh name_fasta
./ipsort_leshy.sh $1
./mitofates_leshy.sh $1
./nommpred_leshy.sh $1
python pts_search.py $1
./targetp.sh $1
./run_multiloc2_with_iprscan.sh $1 animal ${1%.*}'.ML2_animalLO.txt'
./run_multiloc2_with_iprscan.sh $1 fungal ${1%.*}'.ML2_fungalLO.txt'
#pandas required:
python targeting-merger_mito.py -p ${1%.*} -f $1 -d ~/localize