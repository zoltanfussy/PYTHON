# libraries
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
 
# Make a data frame with the GPS of a few cities:
data = pd.DataFrame({
'lat':[-1.383, 15,    6.75,  31.47, 23.564, 42.5,  29.151,  23, 14.449,  40.9,  40.9, 77.833, -77.833, -64.779, 53.5],
'lon':[-89.65, 97, -53.317, -80.42, 58.853,  -68, -85.542, -75, 64.999, 14.15, 14.15, -75.55,     163, -64.058, 4.55],
'name':['CCMP1528', 'CCMP1524', 'CCMP628', 'CCMP2754', 
'CCMP2710', 'CCMP1805', 'CCMP627', 'CCMP629',
'CCMP2000', 'CCMP3104', 'CCMP2496', 'CCMP2192',
'Phaant CCMP1374', 'CCMP1871', "Phaglo PgG"]
})

print(data)
 
#CCMP1528: 1.3833° S  89.65° W
#CCMP1524: 15° N  100° E very approximate
#CCMP628: 6.75° N  53.3167° W
#CCMP2754: 31.47° N  80.42° W
#CCMP2710: 23.5643° N  58.853° E
#CCMP1805: 42.5° N  68° W
#CCMP627: 29.1506° N  85.542° W
#CCMP629: 23° N  75° W  very approximate
#CCMP2000: 14.449° N  64.9997° E
#CCMP3104: 40.9° N  14.15° E
#CCMP2496: 40.9° N  14.15° E
#CCMP2192: 77.8333° N  75.55° W
#CCMP1374: 77.8333° S  163° E 
#CCMP1871: 64.7792° S  64.0575° W 


# A basic map
m=Basemap(llcrnrlon=-170, llcrnrlat=-85,urcrnrlon=189,urcrnrlat=85)
m.drawmapboundary(fill_color='white', linewidth=0) #blueish: #A6CAE0
m.fillcontinents(color='white', alpha=0.7, lake_color='white')
m.drawcoastlines(linewidth=0.8, color="gray")
 
# Add a marker per city of the data frame!
m.plot(data['lon'], data['lat'], label=data["name"], linestyle='none', marker="o", markersize=7, alpha=0.6, c="#A6CAE0", markeredgecolor="black", markeredgewidth=1.5)
for i,j in enumerate(data["name"]):
	plt.text(data['lon'][i]+1, data['lat'][i]+1, j)
plt.savefig("phaeomap.eps")
plt.show()