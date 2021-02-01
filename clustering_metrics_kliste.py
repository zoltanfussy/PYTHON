import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.spatial.distance import cdist,pdist
from scipy.stats import pearsonr
from math import ceil, isnan
import numpy as np
import pandas as pd
import os

os.chdir(wd)

indices = ['egg', 'larvae', 'fednymph', 'female', 'fedfemale']

print("loading data...")
#read pickled data prepared by assign_cluster.py:
read_from_pickle = True
if read_from_pickle:
    X = pd.read_pickle("tick_clustered_expression_matrix.pkl")
    allcolumns = list(X.columns)
    print("identified columns:")
    print(", ".join(allcolumns))
else:
    print("failed")

Xsub = X.loc[:,"egg":"fedfemale"]
treatments = allcolumns[:allcolumns.index("egg")]
print("Found data for clustering methods:")
print(", ".join(treatments))
stagesn = [1, 2, 3, 4, 5]
centroids = {}
reportfile = open("clustering_report_norm_1-Pearson.txt", "w")
errorfile = open("errors_clustering_metrics.txt", "w")
badclusters = set()
plotx = []
ploty = []

allcent = list(Xsub.mean())
print("all points centroid", allcent)
    
for kmean in treatments[:]:   
    cluststats = {}
    maxmodules = int(kmean.split("_")[-1])
    centroids[kmean] = ([], "variability")
    print("Analyzing data for: " + kmean)
    reportfile.write("{} with max modules {}\n".format(kmean, maxmodules))
    cluster_labels_morpheus = np.array(X[kmean], dtype=np.int32)
    for cluster in range(1,maxmodules + 1):
        #determine which contigs belong to this cluster
        contigs = Xsub[X[kmean] == cluster]
        cluststats["{}_{}".format(maxmodules, cluster)] = {"mean": contigs.mean(), 
                                                           "std": contigs.std(), 
                                                           "count": contigs.count()["egg"]}
        # calculate within-cluster variance
        centroids[kmean][0].append(list(contigs.mean()))
    cent = np.array(centroids[kmean][0]) #list of clusters centroids
    
    #calculate minimum distance of each contig to any centroid:    
    D_k = [cdist(Xsub, c, 'euclidean') for c in [cent]]
    #cIdx = [np.argmin(D,axis=1) for D in D_k]
    dist = [np.min(D,axis=1) for D in D_k]
    reportfile.write("distvalues: {}, ".format(len(dist[0])))
    dist = [d[~np.isnan(d)] for d in dist] # ~ >> logical not
    reportfile.write("distvalues without NaN: {}\n".format(len(dist[0])))
    #some distance stats:
    tot_withinss = [sum(d**2) for d in dist]    # Total within-cluster sum of squares
    reportfile.write("total within-cluster distance sqsum:\t{}\n".format(tot_withinss[0]))  
    totss = sum(pdist(Xsub.dropna(how='any'))**2)/Xsub.shape[0]  # The total sum of squares
    reportfile.write("total distance sqsum:\t{}\n".format(totss))
    betweenss = totss - tot_withinss           # The between-cluster sum of squares
    reportfile.write("total explained distance sqsum:\t{}({:.2f}%)\n\n"\
                     .format(betweenss[0], 100 * betweenss[0]/totss))
    
    
    #1-Pearson correlation CALCULATED WITHIN-GROUP/GLOBALLY
    tot_withinps = 0
    tot_ps = 0 # 22741.720695223998 22798.02672928206 #calculated previously from 1-Pearson distance from average centroid
    Xsub_noNaN = Xsub.dropna(how='any')
    c = 0
    for x in range(1, Xsub_noNaN.shape[0] + 1): #go through all valid clusters
        point = np.array(Xsub_noNaN[x-1:x])[0] 
        value = np.min([1 - pearsonr(point, c)[0] for c in cent])
        avgvalue = 1 - pearsonr(point,allcent)[0]
        if not isnan(value):
            tot_withinps += value
            c += 1
            tot_ps += avgvalue
        elif Xsub_noNaN[x-1:x].index[0] not in badclusters:
            errorfile.write("bad cluster {}\n".format(Xsub_noNaN[x-1:x].index[0]))
            badclusters.add(Xsub_noNaN[x-1:x].index[0])

    print("sum of within-cluster variability: {} ({} total)".format(tot_withinps, tot_ps))
    #print(str(c), "valid floats added")
    plotx.append(maxmodules)
    ploty.append((tot_ps - tot_withinps)/tot_ps*100)
    reportfile.write("total within-cluster 1-pearson:\t{}\t{:.2f}% of total explained\n\n"\
                     .format(tot_withinps, (tot_ps - tot_withinps)/tot_ps*100))
    print("Success!")

reportfile.write("table for the elbow plot:\n")
for i in range(len(plotx)):
    reportfile.write("{}\t{}\n".format(plotx[i], ploty[i]))

print("Closing report file, now to create elbow plot...")
reportfile.close()
errorfile.close()


# elbow plot
fig = plt.figure()
ax = fig.add_subplot(111)
#dendrogram clustering is loaded first, remaining are k-means
d = len([x for x in treatments if x.startswith("dendro")])

#dendro values
ax.plot(plotx[:d], ploty[:d], 'bo-', linewidth=0, markersize=6)
#k-means values
ax.plot(plotx[d:], ploty[d:], 'r^--', linewidth=0, markersize=8)

#mark the best cluster number
best = 21 #or 21
#print("best cluster", treatments[best])
ax.plot(plotx[best], ploty[best], marker='o', markersize=12, 
    markeredgewidth=3, markeredgecolor='g', markerfacecolor='None')
ax.set_ylim((60,100)) # 80-100%
plt.grid(True)
plt.xlabel('Number of clusters')
plt.ylabel('Percentage of variance explained (%)')
plt.title('Elbow plot for dendrogram/KMeans 1-Pearson clustering')

plt.savefig("clustering_elbow.pdf") 
print("Analysis done!")

