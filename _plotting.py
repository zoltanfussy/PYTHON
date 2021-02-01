import pandas as pd
brooklyn_one_bed = pd.read_csv('brooklyn-one-bed.csv') #reading tables into pd.array
june = london_data.loc[london_data["month"] == 6]["TemperatureC"] #get temperature for filtered data (june)
#transactions is a pd.dataframe
transactions = transactions.drop(["Unnamed: 0"], axis = 1) #drop a column from data
times = transactions["Transaction Time"].values # data to a separate numpy array
times = transactions["Transaction Time"] #keep column and row info
times = transactions["Transaction Time"].unique() #list of unique items in column
times = transactions[transactions["Transaction Time"] <= whatever] #filter dataframes
times_hist = np.histogram(times, range = (0, 24), bins = 4) #bin data as a histogram
#histogram produces two arrays, first array contains counts per bin, second contains the min (max) values per bin

import numpy as np
#common characteristics of data:
#centrality, modality, spread (range), skew, outliers
average_age = np.average(author_ages) #must be an array
median_age = np.median(author_ages) #must be an array; if the array has even number of items, 
									#it is two values/average of two mid values
#right-skew if average is greater than median, left-skew if opposite, or symmetric distro
example_sum = np.sum(example_array)
example_variance = np.var(example_array) #diff from average squared
example_std = np.std(example_array) #sqrt of variance, i.e. ** 0.5
example_min = np.amin(times) #amax for array max

from scipy import stats
example_mode = stats.mode(example_array) #most frequent value in data
#unimodal has one peak, bimodal/multimodal more peaks; no peaks - uniform
str(example_mode[0][0]) #the value
str(example_mode[1][0]) #its count

data_q3 = np.quantile(example_array, 0.75) #upper quartile value
data_quantiles = np.quantile(example_array, [0.2, 0.4, 0.6, 0.8]) #quintile values; 1.0 not needed
#interquartile range (IQR) ignores the outliers, so you know the range around which your data is centered
from scipy.stats import iqr
interquartile_range = iqr(dataset) #or stats.iqr(dataset)

from matplotlib import pyplot as plt
# Use plt.hist() below
plt.figure(figsize=(15,5)) #define dimensions first!
plt.hist(author_ages, range=(10, 80), alpha = 0.75, bins=14,  edgecolor='black') #tzn 10-15,15-20,20-25,...
#alpha >> transparency
#if more histograms need be plotted, just repeat plt.hist with another data
plt.title("Age of Top 100 Authors at Publication")
plt.axvline(average_age, color='r', linestyle='solid', linewidth=2, label="Mean") #add vertical line
#dotted and dashed lines ale alternatives
plt.xlabel("Age")
plt.ylabel("Count")
# label italic style
ylabel = "$\it{P. antarctica}$ | $\it{P. globosa}$"
plt.ylabel(ylabel, fontsize=20)
plt.set_xlim(minimum, maximum)
#if you have variability in data, use ribbon plot:
plt.fill_between(data, y_lower, y_upper, alpha=0.2)

plt.show()

#if two subplots are needed, then start with 
plt.figure(1) #might not be necessary
plt.figure(figsize=(20,6)) #if you need to define the size of the result
plt.subplot(211)
whatever plot

plt.subplot(212)
whatever plot below that

plt.tight_layout() #for plots less squished

#boxplot is cool because it plots automatically several datasets side-by-side
plt.boxplot([dataset_one, dataset_two], labels = ["Label 1", "Label 2"])
#!!whiskers extend from the box 1.5 times the size of the IQR, outliers are shown as dots
#to quickly fill boxplot datasets from pandas:
states = chest_pain["Provider State"].unique()
datasets = []
for state in states:
  datasets.append(chest_pain[chest_pain['Provider State'] == state][' Average Covered Charges '].values)
plt.boxplot(datasets, labels=states)
plt.savefig("figname.png")
plt.show()

#plotting from pandas dataframe, with x-values added as extra row
#check out on the internet how to plot dataframes with data as columns (pandas default)
rows = list(df.index)
rows.remove("x")
for row in rows:
	#print(df.loc[row])
	if row == "A":
		plt.plot(df.loc["x",:], df.loc[row,:], marker='', color='r', linewidth=2, alpha=0.4)
	else:
		plt.plot(df.loc["x",:], df.loc[row,:], marker='', color=blues.pop(), linewidth=1, alpha=0.7) #pop is random

#functional array of plots from clustering_metrics_kliste.py
numrows = ceil(maxmodules/5)
subplotids = [i for i in range(numrows * 5)]
fig, axs = plt.subplots(nrows=numrows, ncols=5, figsize=(12.6, 0.6*maxmodules), sharex=False,
                        subplot_kw={'xticks': stagesn}) #'yticks': [] to leave empty ax
fig.subplots_adjust(top=0.95, bottom=0.05, left=0.05, right=0.97, hspace=1, wspace=0.3)
fig.suptitle('RSEM across stages for individual clusters ({})'.format(maxmodules))

for i, data in enumerate(cluststats.keys()):
    y = cluststats[data]["mean"]
    err = cluststats[data]["std"]
    axs[i//5,i%5].plot(stagesn, cluststats[data]["mean"], color="darkred", linewidth=3)
    axs[i//5,i%5].axhline(y=0, color='k')
    #axs[i//5,i%5].xticks(np.arange(5), stages)
    suptitle = "{} ({} genes)".format(data, cluststats[data]["count"])
    axs[i//5,i%5].set_title(suptitle)
    axs[i//5,i%5].fill_between(stagesn, y-err, y+err, facecolor="red", alpha=0.2)
    #some other parameters for the shading: edgecolor=, antialiased=True 
    #axs[i//5,i%5].errorbar(stagesn, cluststats[data]["mean"], yerr=cluststats[data]["std"], \
    #                       linewidth=0, elinewidth=1)
    axs[i//5,i%5].set_xticklabels(stages)

    for tick in axs[i//5,i%5].get_xticklabels():
        tick.set_rotation(55)
    subplotids.remove(i)

for i in subplotids:
    fig.delaxes(axs[i//5,i%5])
plt.legend(loc="center right", framealpha = 0, labelspacing=0.2, borderaxespad=0.1, handlelength=0.5, 
                markerscale=2.0, facecolor="grey")
plt.savefig(kmean + "_shade_plots.pdf")                                      
plt.show()
#####
#https://python4astronomers.github.io/plotting/advanced.html
#tricks for multiple axes

#do stuff for pandas subsets!
for i,k in enumerate(df.column.unique())
    groups = df.groupby('column')
    for name, group in groups:
        print(name, group) #group is a df subset!
        axs[i].plot(some_data)


######
#rotate shared xticks:
for data in dataset:
    #...
    xlabels = axs[i//3,i%3].get_xticks() #get_xticklabels if text
    xlabels = [int(x) for x in xlabels]
    axs[3,i%3].set_xticklabels(xlabels, rotation=45, ha="right")



#######
#plot with attached histograms on axes (blobplot-like)
# Create some normally distributed data
mean = [0, 0]
cov = [[1, 1], [1, 2]]
x, y = np.random.multivariate_normal(mean, cov, 3000).T

# Set up the axes with gridspec
fig = plt.figure(figsize=(6, 6))
grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
main_ax = fig.add_subplot(grid[:-1, 1:])
y_hist = fig.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
x_hist = fig.add_subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)

# scatter points on the main axes
main_ax.plot(x, y, 'ok', markersize=3, alpha=0.2)

# histogram on the attached axes
x_hist.hist(x, 40, histtype='stepfilled',
            orientation='vertical', color='gray')
x_hist.invert_yaxis()

y_hist.hist(y, 40, histtype='stepfilled',
            orientation='horizontal', color='gray')
y_hist.invert_xaxis()


#######
#polar plot with offset origin
# Fixing random state for reproducibility
np.random.seed(19680801)

# Compute areas and colors
N = 150
r = 2 * np.random.rand(N)
theta = 2 * np.pi * np.random.rand(N)
area = 200 * r**2
colors = theta

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
c = ax.scatter(theta, r, c=colors, s=area, cmap='hsv', alpha=0.75)

ax.set_rorigin(-2.5)
ax.set_theta_zero_location('W', offset=10)

#################################################################
#https://towardsdatascience.com/a-step-by-step-guide-for-creating-advanced-python-data-visualizations-with-seaborn-matplotlib-1579d6a1a7d0
#create an irregular array of plots => 7+1

freq = ((df_gp.Age.value_counts(normalize = True).reset_index().sort_values(by = 'index').Age)*100).tolist()
number_gp = 7
# freq = the percentage for each age group, and there’re 7 age groups.
def ax_settings(ax, var_name, x_min, x_max):
    ax.set_xlim(x_min,x_max)
    ax.set_yticks([])
    
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.spines['bottom'].set_edgecolor('#444444')
    ax.spines['bottom'].set_linewidth(2)
    
    ax.text(0.02, 0.05, var_name, fontsize=17, fontweight="bold", transform = ax.transAxes) 
    return None
# Manipulate each axes object in the left. Try to tune some parameters and you'll know how each command works.
fig = plt.figure(figsize=(12,7))
gs = gridspec.GridSpec(nrows=number_gp, 
                       ncols=2, 
                       figure=fig, 
                       width_ratios= [3, 1],
                       height_ratios= [1]*number_gp,
                       wspace=0.2, hspace=0.05
                      )
ax = [None]*(number_gp + 1)
features = ['0-17', '18-25', '26-35', '36-45', '46-50', '51-55', '55+']
# Create a figure, partition the figure into 7*2 boxes, set up an ax array to store axes objects, and create a list of age group names.  
for i in range(number_gp):
    ax[i] = fig.add_subplot(gs[i, 0])
    
    ax_settings(ax[i], 'Age: ' + str(features[i]), -1000, 20000)    
    
    sns.kdeplot(data=df_gp[(df_gp.Gender == 'M') & (df_gp.Age == features[i])].Purchase, 
            ax=ax[i], shade=True, color="blue",  bw=300, legend=False)
    sns.kdeplot(data=df_gp[(df_gp.Gender == 'F') & (df_gp.Age == features[i])].Purchase, 
            ax=ax[i], shade=True, color="red",  bw=300, legend=False)
    
    if i < (number_gp - 1): 
        ax[i].set_xticks([])
# this 'for loop' is to create a bunch of axes objects, and link them to GridSpec boxes. Then, we manipulate them with sns.kdeplot() and ax_settings() we just defined.
ax[0].legend(['Male', 'Female'], facecolor='w')
# adding legends on the top axes object     
ax[number_gp] = fig.add_subplot(gs[:, 1])
ax[number_gp].spines['right'].set_visible(False)
ax[number_gp].spines['top'].set_visible(False)
ax[number_gp].barh(features, freq, color='#004c99', height=0.4)
ax[number_gp].set_xlim(0,100)
ax[number_gp].invert_yaxis()
ax[number_gp].text(1.09, -0.04, '(%)', fontsize=10, transform = ax[number_gp].transAxes)   
ax[number_gp].tick_params(axis='y', labelsize = 14)
# manipulate the bar plot on the right. Try to comment out some of the commands to see what they actually do to the bar plot.
plt.show()