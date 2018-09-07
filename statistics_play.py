
# coding: utf-8

# In[48]:

import numpy as np
import pandas as pd

array = {}

with open("test_matrix.tsv") as f:
    for l in f:
        data = l.strip().split("\t")
        array[data[0]] = [float(x) for x in data[1:]]
print("done")


# In[49]:

df = pd.DataFrame(array, index=["egg", "nymph", "fednymph", "female", "fedfemale"])
X = df.T


# In[78]:

print(X.head(10))


# In[57]:

variance = np.var(X, axis=1) #0=x-axis, 1=y-axis
variance = np.var(X.loc[:,:"fedfemale"], axis=1) #limits calculation to columns defined after the comma


# In[58]:

variance.sort_values(ascending=0, inplace=True)
variance.head(10)


# In[77]:

limitvariance = np.var(X[X.nymph > 1], axis=1)
limitvariance.sort_values(ascending=0, inplace=True)
limitvariance.head(10)


# In[ ]:




# In[ ]:



