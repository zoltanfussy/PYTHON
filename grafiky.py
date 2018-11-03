import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol

def get_cmap(n, name='plasma'): #hsv for very divergent data?
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    colormap = plt.cm.get_cmap(name, n)
    rgbcolors = []
    for i in range(colormap.N):
        rgb = colormap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        rgbcolors.append(mcol.rgb2hex(rgb))
    return rgbcolors

#read in some tabular data, or random stuff
y = [np.random.rand(80)*2*171,np.random.rand(81)*2*172,np.random.rand(83)*2*174,np.random.rand(83)*2*172,np.random.rand(163)*2*82] #dosad hodnoty velkosti modulov
#works with lists of numbers too
#y = [list(np.random.rand(80)*2*171),list(np.random.rand(81)*2*172),list(np.random.rand(83)*2*174),list(np.random.rand(83)*2*172),list(np.random.rand(163)*2*82)] #dosad hodnoty velkosti modulov
x = [1,2,3,4,5] #zoznam modulov, musia byt cisla myslim

data = []
for i in y:
    data.append(np.array(i))
fig = plt.figure()
boxprops = dict(linestyle='--', linewidth=2, color='black')
medianprops = dict(linestyle='-.', linewidth=2.5, color='firebrick')
bplot = plt.boxplot(data, notch=True, patch_artist=True, boxprops=boxprops, medianprops=medianprops)

colors = get_cmap(len(x)) #alebo definuj: ['pink', 'lightblue', 'lightgreen']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

for xe, ye in zip(x, y):
    plt.plot([xe] * len(ye), ye, 'o', mfc='none', c='black', zorder=10) #zorder aby boli hore
    
plt.xticks([1, 2, 3, 4, 5])
plt.axes().set_xticklabels(['D.jap', 'R.hum', 'L.lan', 'S.spe', 'YPF'])

plt.show()