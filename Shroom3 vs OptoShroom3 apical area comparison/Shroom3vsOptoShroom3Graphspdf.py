# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 08:50:26 2022

@author: ara
"""

##Document to analyze Shroom3 vs OptoShroom3 apical constriction:

    
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import join


s2hpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hdata'

names2h = listdir(s2hpath)
xyres = 0.2071602

vol2hstim = []
vol2hnonstim = []
relvol2hstim = []
relvol2hnonstim = []
length2hstim = []
length2hnonstim = []
rellength2hstim = []
rellength2hnonstim = []
apical2hstim = []
apical2hnonstim = []
relapical2hstim = []
relapical2hnonstim = []
basal2hstim = []
basal2hnonstim = []
relbasal2hstim = []
relbasal2hnonstim = []
for i in names2h:
    path = join(s2hpath, i)
    temp = np.load(path, allow_pickle = True)
    vol2hstim = vol2hstim + [list(temp[0][a]) for a in range(len(temp[0]))]
    vol2hnonstim = vol2hnonstim + [list(temp[1][a]) for a in range(len(temp[1]))]
    relvol2hstim = relvol2hstim + [list(temp[2][a]) for a in range(len(temp[2]))]
    relvol2hnonstim = relvol2hnonstim + [list(temp[3][a]) for a in range(len(temp[3]))]
    length2hstim = length2hstim + [list(temp[16][a]) for a in range(len(temp[16]))]
    length2hnonstim = length2hnonstim + [list(temp[17][a]) for a in range(len(temp[17]))]
    rellength2hstim = rellength2hstim + [list(temp[18][a]) for a in range(len(temp[18]))]
    rellength2hnonstim = rellength2hnonstim + [list(temp[19][a]) for a in range(len(temp[19]))]
    apical2hstim = apical2hstim + [list(temp[8][a]) for a in range(len(temp[8]))]
    apical2hnonstim = apical2hnonstim + [list(temp[9][a]) for a in range(len(temp[9]))]
    relapical2hstim = relapical2hstim + [list(temp[10][a]) for a in range(len(temp[10]))]
    relapical2hnonstim = relapical2hnonstim + [list(temp[11][a]) for a in range(len(temp[11]))]
    basal2hstim = basal2hstim + [list(temp[12][a]) for a in range(len(temp[12]))]
    basal2hnonstim = basal2hnonstim + [list(temp[13][a]) for a in range(len(temp[13]))]
    relbasal2hstim = relbasal2hstim + [list(temp[14][a]) for a in range(len(temp[14]))]
    relbasal2hnonstim = relbasal2hnonstim + [list(temp[15][a]) for a in range(len(temp[15]))]
    
#Correction: Original segmentation code did output the length in pixels, we need it in microns:

length2hstim = [[a*xyres for a in i] for i in length2hstim]
length2hnonstim = [[a*xyres for a in i] for i in length2hnonstim]


#I decided to make a function to make averages of lists

def getAverage(nestedlist):
    #For 2 layer nested lists of same lenght == repetitions of same measurements of different experiments
    averagelist = []
    standevup = []
    standevlow = []
    for i in range(len(nestedlist[0])):
        temp = []
        for n in range(len(nestedlist)):
            temp.append(nestedlist[n][i])
        averagelist.append(np.mean(temp))
        standevup.append(np.mean(temp) + np.std(temp))
        standevlow.append(np.mean(temp) - np.std(temp))
    return [averagelist, standevup, standevlow]

vol2hstimls = getAverage(vol2hstim)
vol2hnonstimls = getAverage(vol2hnonstim)
relvol2hstimls = getAverage(relvol2hstim)
relvol2hnonstimls = getAverage(relvol2hnonstim)
length2hstimls = getAverage(length2hstim)
length2hnonstimls = getAverage(length2hnonstim)
rellength2hstimls = getAverage(rellength2hstim)
rellength2hnonstimls = getAverage(rellength2hnonstim)
apical2hstimls = getAverage(apical2hstim)
apical2hnonstimls = getAverage(apical2hnonstim)
relapical2hstimls = getAverage(relapical2hstim)
relapical2hnonstimls = getAverage(relapical2hnonstim)
basal2hstimls = getAverage(basal2hstim)
basal2hnonstimls = getAverage(basal2hnonstim)
relbasal2hstimls = getAverage(relbasal2hstim)
relbasal2hnonstimls = getAverage(relbasal2hnonstim)


#Now for OptoShroom3 we need to get the values at timepoint 3 (before stimulation) and timepoint 15 (last point of stimulation)
apical2hstimls[0][15]
apsOptoOff = []
apsOptoOn = []
for i in range(len(apical2hstim)):
    apsOptoOff.append(apical2hstim[i][3])
    apsOptoOn.append(apical2hstim[i][15])

#Same for basal!!
basOptoOff = []
basOptoOn = []
for i in range(len(basal2hstim)):
    basOptoOff.append(basal2hstim[i][3])
    basOptoOn.append(basal2hstim[i][15])

#Let's get apicobasal ratio:
ratioOptoOff = []
for i in range(len(apsOptoOff)):
    ratioOptoOff.append(apsOptoOff[i]/basOptoOff[i])
    

ratioOptoOn = []
for i in range(len(apsOptoOn)):
    ratioOptoOn.append(apsOptoOn[i]/basOptoOn[i])

#Now volumes:
volOptoOff = []
volOptoOn = []
for i in range(len(vol2hstim)):
    volOptoOff.append(vol2hstim[i][3])
    volOptoOn.append(vol2hstim[i][15])

#And lengths:
lenOptoOff = []
lenOptoOn = []
for i in range(len(length2hstim)):
    lenOptoOff.append(length2hstim[i][3])
    lenOptoOn.append(length2hstim[i][15])

#Perfect! Now we need to gather the Shroom3 and GFPCAAX data:

shrpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\Shroom3VolumeMeasurements\ShroomData'
gfppath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\Shroom3VolumeMeasurements\GFPData'
graphpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\Shroom3VolumeMeasurements\Graphs'

namesshr= listdir(shrpath)
namesgfp= listdir(gfppath)

#Shroom3
apsshr = []
basshr = []
lenshr = []
volshr = []
for i in namesshr:
    path = join(shrpath, i)
    temp = np.load(path, allow_pickle = True)
    volshr = volshr + list(temp[0,:])
    apsshr = apsshr + list(temp[1,:])
    basshr = basshr + list(temp[2,:])
    lenshr = lenshr + list(temp[3,:])

#GFP:

apsgfp = []
basgfp = []
lengfp = []
volgfp = []
for i in namesgfp:
    path = join(gfppath, i)
    temp = np.load(path, allow_pickle = True)
    
    volgfp = volgfp + list(temp[0])
    apsgfp = apsgfp + list(temp[1])
    basgfp = basgfp + list(temp[2])
    lengfp = lengfp + list(temp[3])

#Let's quickly get ratios:

ratiogfp = []
for i in range(len(apsgfp)-1):
    ratiogfp.append(apsgfp[i]/basgfp[i]) #Not sure we can trust this ratio!!!

ratioshr = []
for i in range(len(apsshr)):
    ratioshr.append(apsshr[i]/basshr[i])

###Time for plotting!!!

import seaborn as sns
import pandas as pd 



plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams.update({'font.size': 7})
plt.rcParams['pdf.fonttype'] = 42
# Create an array with the colors you want to use
colors = ["#4f81c2", "#dc756d", "#820263", "#60D394"]
# Set your custom color palette
sns.set_palette(sns.color_palette(colors))

#Including ttest between OptoShroom3 and Shroom3:
from scipy.stats import shapiro, ttest_ind
tstat, pval = ttest_ind(apsOptoOn, apsshr)

print(pval)


plt.figure(figsize = (4.3, 4.3))
sample1 = apsOptoOff
sample2 = apsOptoOn
sample3 = apsshr
sample4 = apsgfp
data = pd.DataFrame({'OptoShroom3 \n non-stimulated': pd.Series(sample1), 'OptoShroom3 \n stimulated': pd.Series(sample2), 'Shroom3': pd.Series(sample3), 'WT control': pd.Series(sample4)})
#top = max(sample1+sample2+sample3+sample4)
#bottom = min(sample1+sample2+sample3+sample4)
ax = sns.swarmplot(data=data, size=3)
ax = sns.boxplot(data=data, boxprops={'facecolor':'None'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel("Apical area " + r"$(μm^{2})$")
#ax.set_ylim( bottom+ 0.1*bottom, top + 0.3*top)
ax.text(1.15, 254, "p = " + "%.2e"%pval )
ax.plot([1,2],[250,250], c='k')
plt.savefig(graphpath  + r"\\apicalareaALL.pdf", bbox_inches='tight', transparent=True)
plt.show()


len(apsshr)
len(apsgfp)

#Export source data:
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"
dataexp = pd.DataFrame({'OptoShroom3 non-stimulated': pd.Series(apsOptoOff), 'OptoShroom3 stimulated': pd.Series(apsOptoOn), 'Shroom3': pd.Series(apsshr), 'WT control': pd.Series(apsgfp),})

dataexp.to_csv(filepath + '\\Shroom3apical.csv', index = False)
