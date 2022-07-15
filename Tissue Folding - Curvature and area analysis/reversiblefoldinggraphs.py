# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:19:58 2020

@author: ara
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import date
today = date.today()
d1 = today.strftime("%Y_%m_%d")

from os.path import join
DataOutputFolder = r"C:\Users\ara\Documents\Python Scripts\OptoShroom3 folding matrigel beads\foldingoutput"
FolderToSave = r"C:\Users\ara\Documents\Python Scripts\OptoShroom3 folding matrigel beads\Graphs\foldinggraphs"
from os import listdir

#Make plot font bigger!!
plt.rcParams.update({'font.size': 45})

listofnames = listdir(DataOutputFolder)
print(listofnames)
NumbColonies = len(listofnames)
print(NumbColonies)

listoffiles = []
for i in range(len(listofnames)):
    loadfile = join(DataOutputFolder, listofnames[i])
    temp = np.load(loadfile, allow_pickle=True)
    listoffiles.append(temp)

# Cool! from here we will have all the data to play with! 
print(listoffiles[0][2])


time = 19 
#Area, ax size and curvature averages of main axis, clasiffied by C or Ind
listofareasI = []
listofsizesI = []
listofcurvI = []
listofcurvIm = [] #for minor axis
for i in range(NumbColonies):

    if i ==2 or i == 3 :
        continue
    else:
        print(i)
        print(len(listoffiles[i][5]))
        listofareasI.append(listoffiles[i][5])
        listofsizesI.append(listoffiles[i][2])
        listofcurvI.append(listoffiles[i][6])
        listofcurvIm.append(listoffiles[i][7])


#Cool now let's do the averages and sd:

areasI = []
areasIplussd = []
areasIminussd = []
sizesI = []
curvI = []
curvIplussd = []
curvIminussd = []
curvIm = []
curvImplussd = []
curvImminussd = []
for t in range(time):
    aI =[]
    sI =[]
    cI =[]
    cIm = []
    for i in range(len(listofareasI)):
        aI.append(listofareasI[i][t]/listofareasI[i][0]) #Normalizing data!
        sI.append(listofsizesI[i][t]/listofsizesI[i][0]) #Normalizing data!
        cI.append(listofcurvI[i][t]) #Not for curvature!
        cIm.append(listofcurvIm[i][t]) 
    areasI.append(np.mean(aI))
    areasIplussd.append(np.mean(aI) + np.std(aI))
    areasIminussd.append(np.mean(aI) - np.std(aI))
    sizesI.append(np.mean(sI))
    curvI.append(np.mean(cI))
    curvIplussd.append(np.mean(cI) + np.std(cI))
    curvIminussd.append(np.mean(cI) - np.std(cI))
    curvIm.append(np.mean(cIm))
    curvImplussd.append(np.mean(cIm) + np.std(cIm))
    curvImminussd.append(np.mean(cIm) - np.std(cIm))

###CooooL! Time to make some graphs!
tp = np.linspace(0, time-1, time)
tp = tp.astype(int)
realtime = tp*2

# change font
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

fig,ax = plt.subplots(figsize = (14, 14), dpi=200)
#plt.title('C450V mutant')
ax.plot(realtime,areasI, c="#dc756d", label="Stimulated", linewidth=6.0)
ax.axvspan(0, 24, color = '#78e3d8', alpha = 0.3)
ax.fill_between(realtime, areasIplussd, areasIminussd, color='#dc756d', alpha=.2)
ax.set_xlabel('Time (h)')
ax.set_ylabel('Relative projected colony area')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig(FolderToSave + "/" + d1+  'MeanAreasfig.png', bbox_inches='tight')
plt.show()


#Number of samples used:
len(listofareasI)

#Exporting source data:
import pandas as pd
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"
dataexp = pd.DataFrame(list(zip(realtime, areasI, areasIplussd, areasIminussd)), columns = ["Time (h)", "Average area stimulated", "Average + standard deviation", "Average - standard deviation"])

dataexp.to_csv(filepath + '\\reversiblefolding.csv', index = False)
    

