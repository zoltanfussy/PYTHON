#!/usr/bin/python
from wordcloud import WordCloud, STOPWORDS
from PIL import Image
import os
import numpy as np
from Bio import Entrez, Medline


text = ""

MAX_COUNT = 50
TERM = 'Füssy Zoltán'

Entrez.email = 'zoltan.fussy@gmail.com' # put your mail here
h = Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=TERM)
result = Entrez.read(h)
ids = result['IdList']
h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
records = Medline.parse(h)

for paper in records:
    try:
        if paper['AB']:
            text = text + ' ' + paper['AB']
    except:
        pass

def create_wordcloud(text):
    mask = np.array(Image.open('wordcloud.png'))
    stopwords = set(STOPWORDS)

    wc = WordCloud(background_color="white", mask=mask,
        max_words = 50, stopwords=stopwords)

    wc.generate(text)

    wc.to_file('Zoltán.png')

create_wordcloud(text)