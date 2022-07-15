# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 16:35:12 2022

@author: ara
"""

#Code for the analysis and graphs of previously 3D segmented data from MDCK optoshroom3 iRFPCAAX cells

#First some code to load the data
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import join

s20minpath = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\20mindata'
s2hpath = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hdata'
graphpath = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\finalgraphs\final2'
#This will be the way to load this data
#temp = np.load(DataOut, allow_pickle = True)
names20min = listdir(s20minpath)
names2h = listdir(s2hpath)

xyres = 0.2071602

vol20minstim = []
vol20minnonstim = []
relvol20minstim = []
relvol20minnonstim = []
length20minstim = []
length20minnonstim = []
rellength20minstim = []
rellength20minnonstim = []
apical20minstim = []
apical20minnonstim = []
relapical20minstim = []
relapical20minnonstim = []
basal20minstim = []
basal20minnonstim = []
relbasal20minstim = []
relbasal20minnonstim = []
for i in names20min:
    path = join(s20minpath, i)
    temp = np.load(path, allow_pickle = True)
    vol20minstim = vol20minstim + temp[0]
    vol20minnonstim = vol20minnonstim + temp[1]
    relvol20minstim = relvol20minstim + temp[2]
    relvol20minnonstim = relvol20minnonstim + temp[3]
    length20minstim = length20minstim + temp[16]
    length20minnonstim = length20minnonstim + temp[17]
    rellength20minstim = rellength20minstim + temp[18]
    rellength20minnonstim = rellength20minnonstim + temp[19]
    apical20minstim = apical20minstim + temp[8]
    apical20minnonstim = apical20minnonstim + temp[9]
    relapical20minstim = relapical20minstim + temp[10]
    relapical20minnonstim = relapical20minnonstim + temp[11]
    basal20minstim = basal20minstim + temp[12]
    basal20minnonstim = basal20minnonstim + temp[13]
    relbasal20minstim = relbasal20minstim + temp[14]
    relbasal20minnonstim = relbasal20minnonstim + temp[15]

#Correction: Original segmentation code did output the length in pixels, we need it in microns:

length20minstim = [[a*xyres for a in i] for i in length20minstim]
length20minnonstim = [[a*xyres for a in i] for i in length20minnonstim]




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

#Good! Now make the averages and get the graphs!
vol20minstimls = getAverage(vol20minstim)
vol20minnonstimls = getAverage(vol20minnonstim)
relvol20minstimls = getAverage(relvol20minstim)
relvol20minnonstimls = getAverage(relvol20minnonstim)
length20minstimls = getAverage(length20minstim)
length20minnonstimls = getAverage(length20minnonstim)
rellength20minstimls = getAverage(rellength20minstim)
rellength20minnonstimls = getAverage(rellength20minnonstim)
apical20minstimls = getAverage(apical20minstim)
apical20minnonstimls = getAverage(apical20minnonstim)
relapical20minstimls = getAverage(relapical20minstim)
relapical20minnonstimls = getAverage(relapical20minnonstim)
basal20minstimls = getAverage(basal20minstim)
basal20minnonstimls = getAverage(basal20minnonstim)
relbasal20minstimls = getAverage(relbasal20minstim)
relbasal20minnonstimls = getAverage(relbasal20minnonstim)

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



#OK! Now on to making graphs!
stimstart = 30
stimend20min = 50
stimend2h = 150
timeres = 10 #10 mins/60mins ( in hours)
frames20min = [i*timeres for i in range(len(vol20minstimls[0]))]
frames2h = [i*timeres for i in range(len(vol2hstimls[0]))]


#Let's make some t-test on stimulated cells, for all 3 parameters, and compare with non stimulated cells
#I will do it all with total values: 

#We need to get for each exp condition, for each parameter, and both for stimulated and non stimulated:
#The first value before stimulation and the last  value of the stimulation
firstval = 3
lastval20min = 5 #last timepoint of stimulation
lastval2h = 15

firstvol20minstim = []
lastvol20minstim = []
firstlength20minstim = []
lastlength20minstim = []
firstapical20minstim = []
lastapical20minstim = []
firstbasal20minstim = []
lastbasal20minstim = []
for i in range(len(vol20minstim)):
    firstvol20minstim.append(vol20minstim[i][firstval])
    lastvol20minstim.append(vol20minstim[i][lastval20min])
    firstlength20minstim.append(length20minstim[i][firstval])
    lastlength20minstim.append(length20minstim[i][lastval20min])
    firstapical20minstim.append(apical20minstim[i][firstval])
    lastapical20minstim.append(apical20minstim[i][lastval20min])
    firstbasal20minstim.append(basal20minstim[i][firstval])
    lastbasal20minstim.append(basal20minstim[i][lastval20min])



firstvol20minnonstim = []
lastvol20minnonstim = []
firstlength20minnonstim = []
lastlength20minnonstim = []
firstapical20minnonstim = []
lastapical20minnonstim = []
firstbasal20minnonstim = []
lastbasal20minnonstim = []
for i in range(len(vol20minnonstim)):
    firstvol20minnonstim.append(vol20minnonstim[i][firstval])
    lastvol20minnonstim.append(vol20minnonstim[i][lastval20min])
    firstlength20minnonstim.append(length20minnonstim[i][firstval])
    lastlength20minnonstim.append(length20minnonstim[i][lastval20min])
    firstapical20minnonstim.append(apical20minnonstim[i][firstval])
    lastapical20minnonstim.append(apical20minnonstim[i][lastval20min])
    firstbasal20minnonstim.append(basal20minnonstim[i][firstval])
    lastbasal20minnonstim.append(basal20minnonstim[i][lastval20min])

#And 2h values!

firstvol2hstim = []
lastvol2hstim = []
firstlength2hstim = []
lastlength2hstim = []
firstapical2hstim = []
lastapical2hstim = []
firstbasal2hstim = []
lastbasal2hstim = []
for i in range(len(vol2hstim)):
    firstvol2hstim.append(vol2hstim[i][firstval])
    lastvol2hstim.append(vol2hstim[i][lastval2h])
    firstlength2hstim.append(length2hstim[i][firstval])
    lastlength2hstim.append(length2hstim[i][lastval2h])
    firstapical2hstim.append(apical2hstim[i][firstval])
    lastapical2hstim.append(apical2hstim[i][lastval2h])
    firstbasal2hstim.append(basal2hstim[i][firstval])
    lastbasal2hstim.append(basal2hstim[i][lastval2h])



firstvol2hnonstim = []
lastvol2hnonstim = []
firstlength2hnonstim = []
lastlength2hnonstim = []
firstapical2hnonstim = []
lastapical2hnonstim = []
firstbasal2hnonstim = []
lastbasal2hnonstim = []
for i in range(len(vol2hnonstim)):
    firstvol2hnonstim.append(vol2hnonstim[i][firstval])
    lastvol2hnonstim.append(vol2hnonstim[i][lastval2h])
    firstlength2hnonstim.append(length2hnonstim[i][firstval])
    lastlength2hnonstim.append(length2hnonstim[i][lastval2h])
    firstapical2hnonstim.append(apical2hnonstim[i][firstval])
    lastapical2hnonstim.append(apical2hnonstim[i][lastval2h])
    firstbasal2hnonstim.append(basal2hnonstim[i][firstval])
    lastbasal2hnonstim.append(basal2hnonstim[i][lastval2h])


#Now for C450V controls
C450Vpath = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\C450Vdata'
stimendC450V = 150

namesC450V = listdir(C450Vpath)
volC450Vstim = []
volC450Vnonstim = []
relvolC450Vstim = []
relvolC450Vnonstim = []
lengthC450Vstim = []
lengthC450Vnonstim = []
rellengthC450Vstim = []
rellengthC450Vnonstim = []
apicalC450Vstim = []
apicalC450Vnonstim = []
relapicalC450Vstim = []
relapicalC450Vnonstim = []
basalC450Vstim = []
basalC450Vnonstim = []
relbasalC450Vstim = []
relbasalC450Vnonstim = []
for i in namesC450V:
    path = join(C450Vpath, i)
    temp = np.load(path, allow_pickle = True)
    volC450Vstim = volC450Vstim + [list(temp[0][a]) for a in range(len(temp[0]))]
    volC450Vnonstim = volC450Vnonstim + [list(temp[1][a]) for a in range(len(temp[1]))]
    relvolC450Vstim = relvolC450Vstim + [list(temp[2][a]) for a in range(len(temp[2]))]
    relvolC450Vnonstim = relvolC450Vnonstim + [list(temp[3][a]) for a in range(len(temp[3]))]
    lengthC450Vstim = lengthC450Vstim + [list(temp[16][a]) for a in range(len(temp[16]))]
    lengthC450Vnonstim = lengthC450Vnonstim + [list(temp[17][a]) for a in range(len(temp[17]))]
    rellengthC450Vstim = rellengthC450Vstim + [list(temp[18][a]) for a in range(len(temp[18]))]
    rellengthC450Vnonstim = rellengthC450Vnonstim + [list(temp[19][a]) for a in range(len(temp[19]))]
    apicalC450Vstim = apicalC450Vstim + [list(temp[8][a]) for a in range(len(temp[8]))]
    apicalC450Vnonstim = apicalC450Vnonstim + [list(temp[9][a]) for a in range(len(temp[9]))]
    relapicalC450Vstim = relapicalC450Vstim + [list(temp[10][a]) for a in range(len(temp[10]))]
    relapicalC450Vnonstim = relapicalC450Vnonstim + [list(temp[11][a]) for a in range(len(temp[11]))]
    basalC450Vstim = basalC450Vstim + [list(temp[12][a]) for a in range(len(temp[12]))]
    basalC450Vnonstim = basalC450Vnonstim + [list(temp[13][a]) for a in range(len(temp[13]))]
    relbasalC450Vstim = relbasalC450Vstim + [list(temp[14][a]) for a in range(len(temp[14]))]
    relbasalC450Vnonstim = relbasalC450Vnonstim + [list(temp[15][a]) for a in range(len(temp[15]))]

#Correction: Original segmentation code did output the length in pixels, we need it in microns:

lengthC450Vstim = [[a*xyres for a in i] for i in lengthC450Vstim]
lengthC450Vnonstim = [[a*xyres for a in i] for i in lengthC450Vnonstim]

volC450Vstimls = getAverage(volC450Vstim)
volC450Vnonstimls = getAverage(volC450Vnonstim)
relvolC450Vstimls = getAverage(relvolC450Vstim)
relvolC450Vnonstimls = getAverage(relvolC450Vnonstim)
lengthC450Vstimls = getAverage(lengthC450Vstim)
lengthC450Vnonstimls = getAverage(lengthC450Vnonstim)
rellengthC450Vstimls = getAverage(rellengthC450Vstim)
rellengthC450Vnonstimls = getAverage(rellengthC450Vnonstim)
apicalC450Vstimls = getAverage(apicalC450Vstim)
apicalC450Vnonstimls = getAverage(apicalC450Vnonstim)
relapicalC450Vstimls = getAverage(relapicalC450Vstim)
relapicalC450Vnonstimls = getAverage(relapicalC450Vnonstim)
basalC450Vstimls = getAverage(basalC450Vstim)
basalC450Vnonstimls = getAverage(basalC450Vnonstim)
relbasalC450Vstimls = getAverage(relbasalC450Vstim)
relbasalC450Vnonstimls = getAverage(relbasalC450Vnonstim)

framesC450V = [i*timeres for i in range(len(volC450Vstimls[0]))]

firstval = 3
lastvalC450V = 15

#And C450V values!

firstvolC450Vstim = []
lastvolC450Vstim = []
firstlengthC450Vstim = []
lastlengthC450Vstim = []
firstapicalC450Vstim = []
lastapicalC450Vstim = []
firstbasalC450Vstim = []
lastbasalC450Vstim = []
for i in range(len(volC450Vstim)):
    firstvolC450Vstim.append(volC450Vstim[i][firstval])
    lastvolC450Vstim.append(volC450Vstim[i][lastvalC450V])
    firstlengthC450Vstim.append(lengthC450Vstim[i][firstval])
    lastlengthC450Vstim.append(lengthC450Vstim[i][lastvalC450V])
    firstapicalC450Vstim.append(apicalC450Vstim[i][firstval])
    lastapicalC450Vstim.append(apicalC450Vstim[i][lastvalC450V])
    firstbasalC450Vstim.append(basalC450Vstim[i][firstval])
    lastbasalC450Vstim.append(basalC450Vstim[i][lastvalC450V])



firstvolC450Vnonstim = []
lastvolC450Vnonstim = []
firstlengthC450Vnonstim = []
lastlengthC450Vnonstim = []
firstapicalC450Vnonstim = []
lastapicalC450Vnonstim = []
firstbasalC450Vnonstim = []
lastbasalC450Vnonstim = []
for i in range(len(volC450Vnonstim)):
    firstvolC450Vnonstim.append(volC450Vnonstim[i][firstval])
    lastvolC450Vnonstim.append(volC450Vnonstim[i][lastvalC450V])
    firstlengthC450Vnonstim.append(lengthC450Vnonstim[i][firstval])
    lastlengthC450Vnonstim.append(lengthC450Vnonstim[i][lastvalC450V])
    firstapicalC450Vnonstim.append(apicalC450Vnonstim[i][firstval])
    lastapicalC450Vnonstim.append(apicalC450Vnonstim[i][lastvalC450V])
    firstbasalC450Vnonstim.append(basalC450Vnonstim[i][firstval])
    lastbasalC450Vnonstim.append(basalC450Vnonstim[i][lastvalC450V])

#All together!

#IMAGES FOR FIGURE 2 of paper:
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['xtick.labelsize']=45
plt.rcParams['ytick.labelsize']=35
plt.rcParams.update({'font.size': 45})
    
fig, axes = plt.subplots(2,2, figsize = (24,24), dpi=200)
row = 0
col = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, apical2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, apical2hstimls[1], apical2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, apical2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, apical2hnonstimls[1], apical2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row].title('Cell apical area')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Cell apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,200)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
row = 0
col = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, basal2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, basal2hstimls[1], basal2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, basal2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, basal2hnonstimls[1], basal2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Cell basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,200)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
row = 1
col = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, length2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, length2hstimls[1], length2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, length2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, length2hnonstimls[1], length2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Cell height " + r"$(μm)$")
axes[row, col].set_ylim(0,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
row = 1
col = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, vol2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, vol2hstimls[1], vol2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, vol2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, vol2hnonstimls[1], vol2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Cell volume (fL)")
axes[row, col].set_ylim(0,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.3, 
                    hspace=0.3)
fig.savefig(graphpath + r"fig23dparam.png", bbox_inches='tight')
plt.show()


# Here we can export this data of figure 2 for source data:
import pandas as pd
sourcepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"


#First data for the figure 2e
dataexp = pd.DataFrame(list(zip(frames2h, apical2hstimls[0], apical2hstimls[2], apical2hstimls[1], apical2hnonstimls[0], apical2hnonstimls[2], apical2hnonstimls[1])),
                       columns = ["Time (min)", "Average apical area stim", "Average apical stim + standard deviation", "Average apical stim - standard deviation",
                                 "Average apical area outer", "Average apical outer + standard deviation", "Average apical outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\fig2ApicalAreas.csv', index = False)
#data for the figure 2f
dataexp = pd.DataFrame(list(zip(frames2h, basal2hstimls[0], basal2hstimls[2], basal2hstimls[1], basal2hnonstimls[0], basal2hnonstimls[2], basal2hnonstimls[1])),
                       columns = ["Time (min)", "Average basal area stim", "Average basal stim + standard deviation", "Average basal stim - standard deviation",
                                 "Average basal area outer", "Average basal outer + standard deviation", "Average basal outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\fig2basalAreas.csv', index = False)

#data for the figure 2g
dataexp = pd.DataFrame(list(zip(frames2h, length2hstimls[0], length2hstimls[2], length2hstimls[1], length2hnonstimls[0], length2hnonstimls[2], length2hnonstimls[1])),
                       columns = ["Time (min)", "Average length stim", "Average length stim + standard deviation", "Average length stim - standard deviation",
                                 "Average length outer", "Average length outer + standard deviation", "Average length outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\fig2lengths.csv', index = False)

#data for the figure 2h
dataexp = pd.DataFrame(list(zip(frames2h, vol2hstimls[0], vol2hstimls[2], vol2hstimls[1], vol2hnonstimls[0], vol2hnonstimls[2], vol2hnonstimls[1])),
                       columns = ["Time (min)", "Average volume stim", "Average volume stim + standard deviation", "Average volume stim - standard deviation",
                                 "Average volume outer", "Average volume outer + standard deviation", "Average volume outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\fig2volumes.csv', index = False)

#Next option is the same with smaller limits

fig, axes = plt.subplots(4,2, figsize = (15,40), dpi=200)
col = 0
row = 0
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, apical20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, apical20minstimls[1], apical20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, apical20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, apical20minnonstimls[1], apical20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row, col].title('Cell apical area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(70,200)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, apical2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, apical2hstimls[1], apical2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, apical2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, apical2hnonstimls[1], apical2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row, col].title('Cell apical area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(70,200)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 1
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, basal20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, basal20minstimls[1], basal20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, basal20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, basal20minnonstimls[1], basal20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(70,200)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, basal2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, basal2hstimls[1], basal2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, basal2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, basal2hnonstimls[1], basal2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(70,200)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 2
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, length20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, length20minstimls[1], length20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, length20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, length20minnonstimls[1], length20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Height " + r"$(μm)$")
axes[row, col].set_ylim(7,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 2
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, length2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, length2hstimls[1], length2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, length2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, length2hnonstimls[1], length2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell length " + r"$(μm)$")
axes[row, col].set_ylim(7,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 3
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, vol20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, vol20minstimls[1], vol20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, vol20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, vol20minnonstimls[1], vol20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Volume (fL)")
axes[row, col].set_ylim(700,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 3
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, vol2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, vol2hstimls[1], vol2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, vol2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, vol2hnonstimls[1], vol2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell volume (fL)")
axes[row, col].set_ylim(700,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
fig.savefig(graphpath + r"fig23dparamoffset.png", bbox_inches='tight')
plt.show()





##Now relative:
    
fig, axes = plt.subplots(4,2, figsize = (20,40), dpi=200)
col = 0
row = 0
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, relapical20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, relapical20minstimls[1], relapical20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, relapical20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, relapical20minnonstimls[1], apical20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row, col].title('Cell apical area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Apical area")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relapical2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, relapical2hstimls[1], relapical2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, relapical2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, relapical2hnonstimls[1], relapical2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row, col].title('Cell apical area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 1
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, relbasal20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, relbasal20minstimls[1], relbasal20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, relbasal20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, relbasal20minnonstimls[1], relbasal20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Basal area")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relbasal2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, relbasal2hstimls[1], relbasal2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, relbasal2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, relbasal2hnonstimls[1], relbasal2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 2
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, rellength20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, rellength20minstimls[1], rellength20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, rellength20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, rellength20minnonstimls[1], rellength20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Height")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 2
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, rellength2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, rellength2hstimls[1], rellength2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, rellength2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, rellength2hnonstimls[1], rellength2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell length " + r"$(μm)$")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 3
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, relvol20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, relvol20minstimls[1], relvol20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, relvol20minnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, relvol20minnonstimls[1], relvol20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell relvolume')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Volume")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 3
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relvol2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, relvol2hstimls[1], relvol2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, relvol2hnonstimls[0], c ='#4f81c2', label = "Non Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, relvol2hnonstimls[1], relvol2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell volume (fL)")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
fig.savefig(graphpath + r"fig23dparamrelative.png", bbox_inches='tight')
plt.show()


#Now comparison with C450V and 2h:
fig, axes = plt.subplots(4,3, figsize = (30,40), dpi=200)
col = 0
row = 0
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, apical20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, apical20minstimls[1], apical20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, apical20minnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames20min, apical20minnonstimls[1], apical20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row, col].title('Cell apical area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,300)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, apical2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, apical2hstimls[1], apical2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, apical2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, apical2hnonstimls[1], apical2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row, col].title('Cell apical area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,300)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 2
row = 0
axes[row, col].axvspan(stimstart, stimendC450V, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(framesC450V, apicalC450Vstimls[0], c ='#60D394', label = "Stimulated \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, apicalC450Vstimls[1], apicalC450Vstimls[2], color='#60D394', alpha=.1)
axes[row, col].plot(framesC450V, apicalC450Vnonstimls[0], c ='#820263', label = "Outer area \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, apicalC450Vnonstimls[1], apicalC450Vnonstimls[2], color='#820263', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#axes[row, col].title('Cell apical area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,300)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 1
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, basal20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, basal20minstimls[1], basal20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, basal20minnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames20min, basal20minnonstimls[1], basal20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,300)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, basal2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, basal2hstimls[1], basal2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, basal2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, basal2hnonstimls[1], basal2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,300)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 2
row = 1
axes[row, col].axvspan(stimstart, stimendC450V, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(framesC450V, basalC450Vstimls[0], c ='#60D394', label = "Stimulated \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, basalC450Vstimls[1], basalC450Vstimls[2], color='#60D394', alpha=.1)
axes[row, col].plot(framesC450V, basalC450Vnonstimls[0], c ='#820263', label = "Outer area \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, basalC450Vnonstimls[1], basalC450Vnonstimls[2], color='#820263', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,300)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 2
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, length20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, length20minstimls[1], length20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, length20minnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames20min, length20minnonstimls[1], length20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Height " + r"$(μm)$")
axes[row, col].set_ylim(0,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 2
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, length2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, length2hstimls[1], length2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, length2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, length2hnonstimls[1], length2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Height " + r"$(μm)$")
axes[row, col].set_ylim(0,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 2
row = 2
axes[row, col].axvspan(stimstart, stimendC450V, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(framesC450V, lengthC450Vstimls[0], c ='#60D394', label = "Stimulated \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, lengthC450Vstimls[1], lengthC450Vstimls[2], color='#60D394', alpha=.1)
axes[row, col].plot(framesC450V, lengthC450Vnonstimls[0], c ='#820263', label = "Outer area \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, lengthC450Vnonstimls[1], lengthC450Vnonstimls[2], color='#820263', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 34)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell length " + r"$(μm)$")
axes[row, col].set_ylim(0,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 3
axes[row, col].axvspan(stimstart, stimend20min, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames20min, vol20minstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames20min, vol20minstimls[1], vol20minstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames20min, vol20minnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames20min, vol20minnonstimls[1], vol20minnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Volume (fL)")
axes[row, col].set_ylim(0,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 3
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, vol2hstimls[0], c ='#dc756d', label = "Stimulated", linewidth = 6)
axes[row, col].fill_between(frames2h, vol2hstimls[1], vol2hstimls[2], color='#dc756d', alpha=.1)
axes[row, col].plot(frames2h, vol2hnonstimls[0], c ='#4f81c2', label = "Outer area", linewidth = 6)
axes[row, col].fill_between(frames2h, vol2hnonstimls[1], vol2hnonstimls[2], color='#4f81c2', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Volume (fL)")
axes[row, col].set_ylim(0,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 2
row = 3
axes[row, col].axvspan(stimstart, stimendC450V, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(framesC450V, volC450Vstimls[0], c ='#60D394', label = "Stimulated \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, volC450Vstimls[1], volC450Vstimls[2], color='#60D394', alpha=.1)
axes[row, col].plot(framesC450V, volC450Vnonstimls[0], c ='#820263', label = "Outer area \nC450V mutant", linewidth = 6)
axes[row, col].fill_between(framesC450V, volC450Vnonstimls[1], volC450Vnonstimls[2], color='#820263', alpha=.1)
axes[row, col].legend(loc = "lower right", fontsize = 44)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell volume (fL)")
axes[row, col].set_ylim(0,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
fig.savefig(graphpath + r"suppfigC450Vcomparison.png", bbox_inches='tight')
plt.show()  

### Let's export all of this data that is used in supplementary figure 5:
dataexp = pd.DataFrame(list(zip(frames20min, apical20minstimls[0], apical20minstimls[2], apical20minstimls[1], apical20minnonstimls[0], apical20minnonstimls[2], apical20minnonstimls[1])),
                       columns = ["Time (min)", "Average apical area stim", "Average apical stim + standard deviation", "Average apical stim - standard deviation",
                                 "Average apical area outer", "Average apical outer + standard deviation", "Average apical outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_20minApicalAreas.csv', index = False)
#data for the figure 2f
dataexp = pd.DataFrame(list(zip(frames20min, basal20minstimls[0], basal20minstimls[2], basal20minstimls[1], basal20minnonstimls[0], basal20minnonstimls[2], basal20minnonstimls[1])),
                       columns = ["Time (min)", "Average basal area stim", "Average basal stim + standard deviation", "Average basal stim - standard deviation",
                                 "Average basal area outer", "Average basal outer + standard deviation", "Average basal outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_20minbasalAreas.csv', index = False)

#data for the figure 2g
dataexp = pd.DataFrame(list(zip(frames20min, length20minstimls[0], length20minstimls[2], length20minstimls[1], length20minnonstimls[0], length20minnonstimls[2], length20minnonstimls[1])),
                       columns = ["Time (min)", "Average length stim", "Average length stim + standard deviation", "Average length stim - standard deviation",
                                 "Average length outer", "Average length outer + standard deviation", "Average length outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_20minlengths.csv', index = False)

#data for the figure 20min
dataexp = pd.DataFrame(list(zip(frames20min, vol20minstimls[0], vol20minstimls[2], vol20minstimls[1], vol20minnonstimls[0], vol20minnonstimls[2], vol20minnonstimls[1])),
                       columns = ["Time (min)", "Average volume stim", "Average volume stim + standard deviation", "Average volume stim - standard deviation",
                                 "Average volume outer", "Average volume outer + standard deviation", "Average volume outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_20minvolumes.csv', index = False)


####and for the C450V clon:
dataexp = pd.DataFrame(list(zip(framesC450V, apicalC450Vstimls[0], apicalC450Vstimls[2], apicalC450Vstimls[1], apicalC450Vnonstimls[0], apicalC450Vnonstimls[2], apicalC450Vnonstimls[1])),
                       columns = ["Time (min)", "Average apical area stim", "Average apical stim + standard deviation", "Average apical stim - standard deviation",
                                 "Average apical area outer", "Average apical outer + standard deviation", "Average apical outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_C450VApicalAreas.csv', index = False)
#basal c450v
dataexp = pd.DataFrame(list(zip(framesC450V, basalC450Vstimls[0], basalC450Vstimls[2], basalC450Vstimls[1], basalC450Vnonstimls[0], basalC450Vnonstimls[2], basalC450Vnonstimls[1])),
                       columns = ["Time (min)", "Average basal area stim", "Average basal stim + standard deviation", "Average basal stim - standard deviation",
                                 "Average basal area outer", "Average basal outer + standard deviation", "Average basal outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_C450VbasalAreas.csv', index = False)

#height c450v
dataexp = pd.DataFrame(list(zip(framesC450V, lengthC450Vstimls[0], lengthC450Vstimls[2], lengthC450Vstimls[1], lengthC450Vnonstimls[0], lengthC450Vnonstimls[2], lengthC450Vnonstimls[1])),
                       columns = ["Time (min)", "Average length stim", "Average length stim + standard deviation", "Average length stim - standard deviation",
                                 "Average length outer", "Average length outer + standard deviation", "Average length outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_C450Vlengths.csv', index = False)

#data for the figure C450V
dataexp = pd.DataFrame(list(zip(framesC450V, volC450Vstimls[0], volC450Vstimls[2], volC450Vstimls[1], volC450Vnonstimls[0], volC450Vnonstimls[2], volC450Vnonstimls[1])),
                       columns = ["Time (min)", "Average volume stim", "Average volume stim + standard deviation", "Average volume stim - standard deviation",
                                 "Average volume outer", "Average volume outer + standard deviation", "Average volume outer - standard deviation"])

dataexp.to_csv(sourcepath + '\\supfig5_C450Vvolumes.csv', index = False)



#Finally comparation of differences for statistical significance:

#VOLUMES:
vol2hnonstimdiffs = [lastvol2hnonstim[i] - firstvol2hnonstim[i] for i in range(len(firstvol2hnonstim))]
vol2hstimdiffs = [lastvol2hstim[i] - firstvol2hstim[i] for i in range(len(firstvol2hstim))]
volC450Vnonstimdiffs = [lastvolC450Vnonstim[i] - firstvolC450Vnonstim[i] for i in range(len(firstvolC450Vnonstim))]
volC450Vstimdiffs = [lastvolC450Vstim[i] - firstvolC450Vstim[i] for i in range(len(firstvolC450Vstim))]

#LENGTHS
length2hnonstimdiffs = [lastlength2hnonstim[i] - firstlength2hnonstim[i] for i in range(len(firstlength2hnonstim))]
length2hstimdiffs = [lastlength2hstim[i] - firstlength2hstim[i] for i in range(len(firstlength2hstim))]
lengthC450Vnonstimdiffs = [lastlengthC450Vnonstim[i] - firstlengthC450Vnonstim[i] for i in range(len(firstlengthC450Vnonstim))]
lengthC450Vstimdiffs = [lastlengthC450Vstim[i] - firstlengthC450Vstim[i] for i in range(len(firstlengthC450Vstim))]

#APICALS
apical2hnonstimdiffs = [lastapical2hnonstim[i] - firstapical2hnonstim[i] for i in range(len(firstapical2hnonstim))]
apical2hstimdiffs = [lastapical2hstim[i] - firstapical2hstim[i] for i in range(len(firstapical2hstim))]
apicalC450Vnonstimdiffs = [lastapicalC450Vnonstim[i] - firstapicalC450Vnonstim[i] for i in range(len(firstapicalC450Vnonstim))]
apicalC450Vstimdiffs = [lastapicalC450Vstim[i] - firstapicalC450Vstim[i] for i in range(len(firstapicalC450Vstim))]

#BASALS
basal2hnonstimdiffs = [lastbasal2hnonstim[i] - firstbasal2hnonstim[i] for i in range(len(firstbasal2hnonstim))]
basal2hstimdiffs = [lastbasal2hstim[i] - firstbasal2hstim[i] for i in range(len(firstbasal2hstim))]
basalC450Vnonstimdiffs = [lastbasalC450Vnonstim[i] - firstbasalC450Vnonstim[i] for i in range(len(firstbasalC450Vnonstim))]
basalC450Vstimdiffs = [lastbasalC450Vstim[i] - firstbasalC450Vstim[i] for i in range(len(firstbasalC450Vstim))]

#Plotting and measuring significances: #WILL HAVE TO ADD PROPER CONTROLS
from scipy.stats import ttest_ind
def twoSampleTComp(sample1, sample2, title, ylab = "y label here" , val1 = "OptoShroom3", val2 = "OptoShroom3 C450V"):
    #Comparing 4 independent samples with t-test and then plotting. 
    tstat, ps1tos2 = ttest_ind(sample1, sample2)
    
    top = max(sample1+sample2) #top value in the graph, to put text over it
    
    plt.figure(figsize = (15,15))
    plt.title(title)
    plt.boxplot([sample1, sample2])
    plt.ylabel(ylab)
    if ps1tos2 > 0.05:
        plt.text(1.5, top+0.05*top, "n.s." )
    else:
        plt.text(1.3, top+0.05*top, "p = " + "%.2e"%ps1tos2 )
    plt.plot([1,2],[top+0.04*top,top+0.04*top], c='k')
    plt.xticks([1,2], [val1,val2])
    plt.savefig(graphpath + title + r".png", bbox_inches='tight')
    plt.show()  
    return ps1tos2

papicals = twoSampleTComp(apical2hstimdiffs, apicalC450Vstimdiffs, title = "apical_surface_differences_2h_stimulation", ylab = "area difference (um**2)")
pvols = twoSampleTComp(vol2hstimdiffs, volC450Vstimdiffs, title = "cell_volume_differences_2h_stimulation", ylab = "volume difference (um**3)")
plens = twoSampleTComp(length2hstimdiffs, lengthC450Vstimdiffs, title = "cell_length_differences_2h_stimulation", ylab = "length difference (um)")
pbasals = twoSampleTComp(basal2hstimdiffs, apicalC450Vstimdiffs, title = "basal_surface_differences_2h_stimulation", ylab = "area difference (um**2)")


#Make the graph:

import seaborn as sns
import pandas as pd 

plt.rcParams.update({'font.size': 16})
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
# Create an array with the colors you want to use
colors = ["#60D394", "#dc756d"]
# Set your custom color palette
sns.set_palette(sns.color_palette(colors))

plt.figure(figsize = (5, 6), dpi=200)
sample1 = apicalC450Vstimdiffs
sample2 = apical2hstimdiffs 
data = pd.DataFrame({'OptoShroom3 C450V': pd.Series(sample1), 'OptoShroom3': pd.Series(sample2)})
top = max(sample1+sample2)
bottom = min(sample1 + sample2)
ax = sns.swarmplot(data=data, size=5)
ax = sns.boxplot(data=data, boxprops={'facecolor':'None'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel("Apical area difference " + r"$(μm^{2})$")
ax.set_ylim( bottom+ 0.1*bottom, top + 0.3*top)
ax.text(0.2, top+0.08*top, "p = " + "%.2e"%papicals )
ax.plot([0,1],[top+0.04*top,top+0.04*top], c='k')
plt.savefig(graphpath  + r"apicalareaboxes.png", bbox_inches='tight')
plt.show()


plt.figure(figsize = (5, 6), dpi=200)
sample1 = basalC450Vstimdiffs
sample2 = basal2hstimdiffs 
data = pd.DataFrame({'OptoShroom3 C450V': pd.Series(sample1), 'OptoShroom3': pd.Series(sample2)})
top = max(sample1+sample2)
bottom = min(sample1 + sample2)
ax = sns.swarmplot(data=data, size=5)
ax = sns.boxplot(data=data, boxprops={'facecolor':'None'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel("Basal area difference " + r"$(μm^{2})$")
ax.set_ylim( bottom+ 0.1*bottom, top + 0.3*top)
ax.text(0.1, top+0.06*top, "p = " + "%.2e"%pbasals )
ax.plot([0,1],[top+0.04*top,top+0.04*top], c='k')
plt.savefig(graphpath  + r"basalareaboxes.png", bbox_inches='tight')
plt.show()


plt.figure(figsize = (5, 6), dpi=200)
sample1 = lengthC450Vstimdiffs
sample2 = length2hstimdiffs 
data = pd.DataFrame({'OptoShroom3 C450V': pd.Series(sample1), 'OptoShroom3': pd.Series(sample2)})
top = max(sample1+sample2)
bottom = min(sample1 + sample2)
ax = sns.swarmplot(data=data, size=5)
ax = sns.boxplot(data=data, boxprops={'facecolor':'None'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel("Height difference " + r"$(μm)$")
ax.set_ylim( bottom+ 0.1*bottom, top + 0.3*top)
ax.text(0.1, top+0.06*top, "p = " + "%.2e"%plens )
ax.plot([0,1],[top+0.04*top,top+0.04*top], c='k')
plt.savefig(graphpath  + r"heightboxes.png", bbox_inches='tight')
plt.show()

plt.figure(figsize = (5, 6), dpi=200)
sample1 = volC450Vstimdiffs
sample2 = vol2hstimdiffs 
data = pd.DataFrame({'OptoShroom3 C450V': pd.Series(sample1), 'OptoShroom3': pd.Series(sample2)})
top = max(sample1+sample2)
bottom = min(sample1 + sample2)
ax = sns.swarmplot(data=data, size=5)
ax = sns.boxplot(data=data, boxprops={'facecolor':'None'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel("Volume difference (fL)")
ax.set_ylim( bottom+ 0.1*bottom, top + 0.3*top)
ax.text(0.1, top+0.06*top, "p = " + "%.2e"%pvols )
ax.plot([0,1],[top+0.04*top,top+0.04*top], c='k')
plt.savefig(graphpath  + r"volumeboxes.png", bbox_inches='tight')
plt.show()

##Let's export this for source data:
dataexp = pd.DataFrame({"Apical area difference OptoShroom3 stim": pd.Series(apical2hstimdiffs), "Apical area difference OptoShroom3 C450V stim": pd.Series(apicalC450Vstimdiffs)})
dataexp.to_csv(sourcepath + '\\supfig5_apicaldifferences.csv', index = False)
#basal
dataexp = pd.DataFrame({"Basal area difference OptoShroom3 stim": pd.Series(basal2hstimdiffs), "Basal area difference OptoShroom3 C450V stim": pd.Series(basalC450Vstimdiffs)})
dataexp.to_csv(sourcepath + '\\supfig5_basaldifferences.csv', index = False)
#length
dataexp = pd.DataFrame({"Height difference OptoShroom3 stim": pd.Series(length2hstimdiffs), " Height difference OptoShroom3 C450V stim": pd.Series(lengthC450Vstimdiffs)})
dataexp.to_csv(sourcepath + '\\supfig5_lengthdifferences.csv', index = False)
#volume 
dataexp = pd.DataFrame({"Volume difference OptoShroom3 stim": pd.Series(vol2hstimdiffs), "Volume difference OptoShroom3 C450V stim": pd.Series(volC450Vstimdiffs)})
dataexp.to_csv(sourcepath + '\\supfig5_volumedifferences.csv', index = False)


#Predicitons:
    #Predict volume:

def predictVolume(apicals, basals, lengths, iterations = 10000):
    #all apicals, basals and lengths must have the same number of timepoints
    volumes = []
    for t in range(len(apicals)):
        print(t)
        vol = lengths[t]*(apicals[t]+basals[t])/2
        volumes.append(vol)
    return volumes

vols = predictVolume(apical2hstimls[0],basal2hstimls[0], length2hstimls[0])

plt.figure()
plt.plot(frames2h, vol2hstimls[0], c= 'blue', label = 'measurement')
plt.plot(frames2h, vols, c = 'red', label = 'prediction')
plt.xlabel("Time (min)")
plt.ylabel("Cell volume (femtoliters)")
plt.ylim(0,1600)
plt.legend()
plt.axvspan(stimstart, stimend2h, color="blue", alpha = 0.1)
plt.show()
 
#Same thing for prediction of length
def predictLength(apicals, basals, volumes, iterations = 10000):
    #all apicals, basals and lengths must have the same number of timepoints
    lengths = []
    for t in range(len(apicals)):
        print(t)
        length = 2*volumes[t]/(apicals[t]+basals[t]) #The length had been divided in i iterations!
        lengths.append(length)
    return lengths

lens = predictLength(apical2hstimls[0],basal2hstimls[0], vol2hstimls[0])

#Amazing! I want to complement with the measurements of apical and basal:

def predictApical(basals, volumes, lengths, iterations = 10000):
    #all apicals, basals and lengths must have the same number of timepoints
    apicals = []
    for t in range(len(basals)):
        print(t)
        apical = 2*volumes[t]/lengths[t] - basals[t]
        apicals.append(apical)
    return apicals

aps = predictApical(basal2hstimls[0], vol2hstimls[0], length2hstimls[0])

#Should also work for basal!
bas = predictApical(apical2hstimls[0], vol2hstimls[0], length2hstimls[0])


#Proportionally?
relvols = vols/vols[0]
rellens = lens/lens[0]
relaps = aps/aps[0]
relbas = bas/bas[0]

#Also for control:
C450Vvols = predictVolume(apicalC450Vstimls[0],basalC450Vstimls[0], lengthC450Vstimls[0])
C450Vlens = predictLength(apicalC450Vstimls[0],basalC450Vstimls[0], volC450Vstimls[0])
C450Vaps = predictApical(basalC450Vstimls[0], volC450Vstimls[0], lengthC450Vstimls[0])

#Should also work for basal!
C450Vbas = predictApical(apicalC450Vstimls[0], volC450Vstimls[0], lengthC450Vstimls[0])

#Proportionally?
relC450Vvols = C450Vvols/C450Vvols[0]
relC450Vlens = C450Vlens/C450Vlens[0]
relC450Vaps = C450Vaps/C450Vaps[0]
relC450Vbas = C450Vbas/C450Vbas[0]


#for the figures:
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['xtick.labelsize']=45
plt.rcParams['ytick.labelsize']=35
plt.rcParams.update({'font.size': 45})
    
#Trial to represent stimulations truely:
stimpoints =[]
for i in range(12):
    stimpoints = stimpoints + [stimstart + i*10 + 1, stimstart + i*10 + 10]
#Not visible change so wont use it
#for i in range(12):
    #axes.axvspan(stimpoints[i*2],stimpoints[2*i+1], color='#a4c3f6', alpha = 0.2)
#Graph1, apical, and lateral predictions:
fig, axes = plt.subplots(figsize = (10,10), dpi=200)
axes.axvspan(stimstart,stimend2h, color='#a4c3f6', alpha = 0.2)
axes.plot(frames2h, length2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes.plot(frames2h, lens, c ='#820263', label = "Predicted", linewidth = 6)
axes.legend(loc = "lower right", fontsize = 44)
#plt.title('Cell height (apical-basal)')
axes.set_xlabel("Time (min)")
axes.set_ylabel("Cell height " + r"$(μm)$")
axes.set_ylim(0,12)
axes.spines['right'].set_visible(False)
axes.spines['top'].set_visible(False)
fig.savefig(graphpath + r"fig2predictors.png", bbox_inches='tight')
plt.show()

#Let's export this data that appears on figure 2j:
dataexp = pd.DataFrame(list(zip(frames2h, length2hstimls[0], lens)),
                       columns = ["Time (min)", "Average length stim", "length prediction"])

dataexp.to_csv(sourcepath + '\\fig2lengthprediction.csv', index = False)

#Graph2 stimulated and C450V predictions
fig, axes = plt.subplots(4,2, figsize = (20,40), dpi=200)
col = 0
row = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, apical2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, aps, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
axes[row, col].set_ylabel("Apical area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,260)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, apicalC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, C450Vaps, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
axes[row, col].set_ylim(0,260)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, basal2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, bas, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,260)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, basalC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, C450Vbas, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,260)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 2
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, length2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, lens, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Height " + r"$(μm)$")
axes[row, col].set_ylim(0,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 2
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, lengthC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, C450Vlens, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell length " + r"$(μm)$")
axes[row, col].set_ylim(0,12)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 3
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, vol2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, vols, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Volume (fL)")
axes[row, col].set_ylim(0,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 3
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, volC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, C450Vvols, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell volume (fL)")
axes[row, col].set_ylim(0,1600)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
fig.savefig(graphpath + r"fig23dpredictors.png", bbox_inches='tight')
plt.show()

#Let's export all absolute values:
dataexp = pd.DataFrame(list(zip(frames2h, apical2hstimls[0], aps, basal2hstimls[0], bas, length2hstimls[0], lens, vol2hstimls[0], vols)),
                       columns = ["Time (min)", "Average apical area stim", "apical area prediction", "Average basal area stim", "basal area prediction", "Average height stim", "height prediction", "Average volume stim", "volume prediction"])

dataexp.to_csv(sourcepath + '\\suppfig7predictionAbsolutesOptoShr.csv', index = False)

dataexp = pd.DataFrame(list(zip(framesC450V, apicalC450Vstimls[0], C450Vaps, basalC450Vstimls[0], C450Vbas, lengthC450Vstimls[0], C450Vlens, volC450Vstimls[0], C450Vvols)),
                       columns = ["Time (min)", "Average apical area C450V stim", "apical area C450V prediction", "Average basal area C450V stim", "basal area C450V prediction", "Average height C450V stim", "height C450V prediction", "Average volume C450V stim", "volume C450V prediction"])

dataexp.to_csv(sourcepath + '\\suppfig7predictionAbsolutesC450V.csv', index = False)

#Graph3 stimulated and C450V predictions relative
fig, axes = plt.subplots(4,2, figsize = (20,40), dpi=200)
col = 0
row = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relapical2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, relaps, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
axes[row, col].set_ylabel("Relative apical area")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 0
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relapicalC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, relC450Vaps, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relbasal2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, relbas, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Relative basal area")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 1
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relbasalC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, relC450Vbas, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell basal area')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell basal area " + r"$(μm^{2})$")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 2
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, rellength2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, rellens, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel("Relative height")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 2
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, rellengthC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, relC450Vlens, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell height (apical-basal)')
#axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell length " + r"$(μm)$")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 0
row = 3
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relvol2hstimls[0], c ='#dc756d', label = "Measured", linewidth = 6)
axes[row, col].plot(frames2h, relvols, c ='#4f81c2', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
axes[row, col].set_ylabel(" Relative volume")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
col = 1
row = 3
axes[row, col].axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
axes[row, col].plot(frames2h, relvolC450Vstimls[0], c ='#60D394', label = "Measured C450V \nmutant", linewidth = 6)
axes[row, col].plot(frames2h, relC450Vvols, c ='#820263', label = "Predicted", linewidth = 6)
axes[row, col].legend(loc = "lower right", fontsize = 38)
#plt.title('Cell volume')
axes[row, col].set_xlabel("Time (min)")
#axes[row, col].set_ylabel("Cell volume (fL)")
axes[row, col].set_ylim(0,1.3)
axes[row, col].spines['right'].set_visible(False)
axes[row, col].spines['top'].set_visible(False)
fig.savefig(graphpath + r"fig23dpredictorsrel.png", bbox_inches='tight')
plt.show()


#And let's export the relative:
dataexp = pd.DataFrame(list(zip(frames2h, relapical2hstimls[0], relaps, relbasal2hstimls[0], relbas, rellength2hstimls[0], rellens, relvol2hstimls[0], relvols)),
                       columns = ["Time (min)", "Average apical area stim", "apical area prediction", "Average basal area stim", "basal area prediction", "Average height stim", "height prediction", "Average volume stim", "volume prediction"])

dataexp.to_csv(sourcepath + '\\suppfig7predictionRelativeOptoShr.csv', index = False)

dataexp = pd.DataFrame(list(zip(framesC450V, relapicalC450Vstimls[0], relC450Vaps, relbasalC450Vstimls[0], relC450Vbas, rellengthC450Vstimls[0], relC450Vlens, relvolC450Vstimls[0], relC450Vvols)),
                       columns = ["Time (min)", "Average apical area C450V stim", "apical area C450V prediction", "Average basal area C450V stim", "basal area C450V prediction", "Average height C450V stim", "height C450V prediction", "Average volume C450V stim", "volume C450V prediction"])

dataexp.to_csv(sourcepath + '\\suppfig7predictionRelativeC450V.csv', index = False)



#Getting number of cells for figure foots:
len(names2h)
len(vol2hstim)
len(vol2hnonstim)

len(names20min)
len(vol20minstim)
len(vol20minnonstim)

len(namesC450V)
len(volC450Vstim)
len(volC450Vnonstim)
