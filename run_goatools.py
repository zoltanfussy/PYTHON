import urllib.request
import os

home = "/Users/zoliq/ownCloud/"
#home = "/Volumes/zoliq data/OwnCloud/"
wd = home + "progs/PYTHON-OTHER/goatools-master/"
os.chdir(wd)
#LOAD THE GO-BASIC file
"""
url = "http://geneontology.org/ontology/go-basic.obo"
fileName = 'go-basic.obo'
urllib.request.urlretrieve(url, fileName)
print("GO-basic file downloaded.")
"""
"python scripts/find_enrichment.py data/study data/population data/association \
--sections=goatools.test_data.sections.data2018_07_find_enrichment \
--pval_field=fdr_bh \
--outfile=goea_fdr_bh.tsv"

#bonferroni can be used too for FDR

"""
for i in kliste/cl*.txt; do j="${i#kliste\/}" && python scripts/find_enrichment.py $i kliste/population.txt kliste/association.tsv \
--sections=goatools.test_data.sections.data2018_07_find_enrichment \
--pval_field=fdr_bh --outfile=$j; done
"""