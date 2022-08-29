# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 20:54:56 2022

@author: ara
"""

#File for making graphs after the analysis of MyosinII cytoskeletal  stainings

import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import join

onpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\cytoskeleton analysis\MyoIIbONdata'
offpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\cytoskeleton analysis\MyoIIbOFFdata'
graphpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\cytoskeleton analysis\graphs'

#I saved files as [junctionaverage, cytoplasmaverage, ratio, junctionall, cytoplasmall]

nameson = listdir(onpath)

onjuncts = []
oncytos = []
onratios = []
for i in nameson:
    path = join(onpath, i)
    temp = np.load(path, allow_pickle = True)
    onjuncts.append(temp[0])
    oncytos.append(temp[1])
    onratios.append(temp[2])

namesoff = listdir(offpath)
offjuncts = []
offcytos = []
offratios = []
for i in namesoff:
    path = join(offpath, i)
    temp = np.load(path, allow_pickle = True)
    offjuncts.append(temp[0])
    offcytos.append(temp[1])
    offratios.append(temp[2])



import seaborn as sns
import pandas as pd 


from scipy.stats import shapiro, ttest_ind
tstat, pval = ttest_ind(offratios, onratios)

print(pval)

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams.update({'font.size': 7})
plt.rcParams['xtick.labelsize']=6
plt.rcParams['ytick.labelsize']=6
plt.rcParams['pdf.fonttype'] = 42
# Create an array with the colors you want to use
colors = ["#4f81c2", "#dc756d", "#820263", "#60D394"]
# Set your custom color palette
sns.set_palette(sns.color_palette(colors))

plt.figure(figsize = (2, 3.6))
data = pd.DataFrame({'OptoShroom3 \n non-stimulated': pd.Series(offratios), 'OptoShroom3 \n stimulated': pd.Series(onratios)})
top = max(offratios+onratios)
bottom = min(offratios+onratios)
ax = sns.swarmplot(data=data, size=4)
ax = sns.boxplot(data=data, boxprops={'facecolor':'None'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel("Myosin IIb apical signal ratio (junction/cytoplasm)")
ax.set_ylim( bottom- 0.1*bottom, top + 0.1*top)
ax.text(0.2, top+0.05*top, "p = " + "%.2e"%pval )
ax.plot([0,1],[top+0.04*top,top+0.04*top], c='k')
plt.savefig(graphpath  + r"\\apicalsignalmyoIIb.pdf", bbox_inches='tight', transparent=True)
plt.show()

#Export source data:
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"
dataexp = pd.DataFrame(list(zip(offratios, onratios)), columns = ["Non-stimulated", "Stimulated"])

dataexp.to_csv(filepath + '\\cytoskeletal analysis.csv', index = False)