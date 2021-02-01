cd /Users/morpholino/OwnCloud/progs/PYTHON-OTHER/goatools-master

for i in dendrogram_22/cl*txt; do
	j="${i#dendrogram_22\/}"
	python scripts/find_enrichment.py $i dendrogram_22/population.txt dendrogram_22/association.tsv \
	--sections=goatools.test_data.sections.data2018_07_find_enrichment \
	--pval_field=fdr_bh --outfile=$j
done