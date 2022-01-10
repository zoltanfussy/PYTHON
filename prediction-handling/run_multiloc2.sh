#run as: ~/run_multiloc2.sh some.fasta animal some.ML2_animalLO.txt
python /opt/multiloc2/git/MultiLoc2/MultiLoc2/src/multiloc2_prediction.py -fasta=$1 -origin=$2 -result=$3 -predictor=LowRes
