#!/usr/bin/python
import numpy as np
#import pandas as pd
#from os import path
from Bio import Entrez, Medline
from PIL import Image
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator

import matplotlib.pyplot as plt
#never name pythons scripts the same as the modules imported!

def blue_color_func(word, font_size, position, orientation, random_state=None,
                    **kwargs):
    #hsl = hue/saturation/lightness
    #hue in rgb degrees, with red being 0, green 120 and blue 240
    #saturation from  0=gray to 100=full colour 
    #lightness from 0=black to 100=white
    return "hsl(240, 50%%, %d%%)" % random.randint(40, 80)

def create_wordcloud(text):
    mask = np.array(Image.open('ceratium.png'))
    stopwords = set(STOPWORDS)

    wc = WordCloud(background_color="white", 
        mask=mask, max_words = 200, max_font_size=50, 
        stopwords=stopwords, 
        collocations=False,
        #colormap="Oranges_r", #in case image_colors suck
        include_numbers=False,
        repeat=False, relative_scaling=0, #r.s. between 0 and 1
        prefer_horizontal=0.7).generate(text)

    image_colors = ImageColorGenerator(mask) #may give an error, but is handy
    #image_colors = blue_color_func #another way to get a colormap
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
#refresh sometimes...
MAX_COUNT = 50
TERM = 'Füssy Zoltán'

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
            text = text + '\n' + paper['AB']
    except:
        pass
print("Imported from NCBI: {} words.".format(len(text.split())))
with open("abstracts_NCBI_automatic.txt", "w") as result:
    result.write(text)
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