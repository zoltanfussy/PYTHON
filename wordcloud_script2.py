#!/usr/bin/python
import numpy as np
#import pandas as pd
#from os import path
from Bio import Entrez, Medline
from PIL import Image
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator

import matplotlib.pyplot as plt
#never name pythons scripts the same as the modules imported!

def create_wordcloud(text):
    mask = np.array(Image.open('ceratium.png'))
    stopwords = set(STOPWORDS)

    wc = WordCloud(background_color="white", 
        mask=mask, max_words = 200, stopwords=stopwords, max_font_size=50,
        prefer_horizontal=0.6).generate(text)

    image_colors = ImageColorGenerator(mask) #may give an error, but is handy
    plt.figure(figsize=[8,6])
    plt.imshow(wc.recolor(color_func=image_colors), interpolation='bilinear')
    plt.axis("off")
    plt.savefig('Zoli_cer.png')
    plt.show()

#################
###    MAIN   ###
#################

text = ""
"""
MAX_COUNT = 50
TERM = 'Věchtová Pavlína'

Entrez.email = 'zoltan.fussy@gmail.com' # put your mail here
h = Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=TERM)
result = Entrez.read(h)
ids = result['IdList']
h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
records = Medline.parse(h)

print("looking up papers at NCBI...")
for paper in records:
    print(paper['TI'])
    try:
        if paper['AB']: #'TI' for titles, 'AB' for abstracts
            text = text + ' ' + paper['AB']
    except:
        pass
print("Imported from NCBI: {} words.".format(len(text.split())))
"""
NCBIabstracts = open("abstracts_NCBI.txt")
for l in NCBIabstracts:
    text = text + ' ' + l
print("Imported from NCBI: {} words.".format(len(text.split())))
#i had problems with non-ascii characters again
import io

notincluded = io.open("abstracts.txt", "r", encoding="utf-8")
for l in notincluded:
    #print(l[:10])
    text = text + ' ' + l
print("Including offline abstracts: {} words.".format(len(text.split())))

print("abstracts read-in finished, now to drawing...")
create_wordcloud(text)