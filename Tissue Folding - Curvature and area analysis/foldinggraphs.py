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
DataOutputFolder = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\OptoShroom Folding\Cured Data m3"
FolderToSave = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\OptoShroom Folding\Graphs m3"
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

time = 13 
#Area, ax size and curvature averages of main axis, clasiffied by C or Ind
listofareasC = []
listofareasI = []
listofsizesC = []
listofsizesI = []
listofcurvC = []
listofcurvI = []
listofcurvCm = [] # for minor axis
listofcurvIm = [] #for minor axis
for i in range(NumbColonies):
    if listoffiles[i][0] == "C":
        listofareasC.append(listoffiles[i][5])
        listofsizesC.append(listoffiles[i][2])
        listofcurvC.append(listoffiles[i][6])
        listofcurvCm.append(listoffiles[i][7])
    else:
        listofareasI.append(listoffiles[i][5])
        listofsizesI.append(listoffiles[i][2])
        listofcurvI.append(listoffiles[i][6])
        listofcurvIm.append(listoffiles[i][7])

#Cool now let's do the averages and sd:
areasC = []
areasCplussd = []
areasCminussd = []
areasI = []
areasIplussd = []
areasIminussd = []
sizesC = []
sizesI = []
curvC = []
curvCplussd = []
curvCminussd = []
curvCm = []
curvCmplussd = []
curvCmminussd = []
curvI = []
curvIplussd = []
curvIminussd = []
curvIm = []
curvImplussd = []
curvImminussd = []
for t in range(time):
    aC =[]
    aI =[]
    sC =[]
    sI =[]
    cC =[]
    cI =[]
    cCm = []
    cIm = []
    for i in range(len(listofareasC)):
        aC.append(listofareasC[i][t]/listofareasC[i][0]) #Normalizing data!
        sC.append(listofsizesC[i][t]/listofsizesC[i][0]) #Normalizing data!
        cC.append(listofcurvC[i][t]) #Not for curvature!
        cCm.append(listofcurvCm[i][t])
        aI.append(listofareasI[i][t]/listofareasI[i][0]) #Normalizing data!
        sI.append(listofsizesI[i][t]/listofsizesI[i][0]) #Normalizing data!
        cI.append(listofcurvI[i][t]) #Not for curvature!
        cIm.append(listofcurvIm[i][t]) 
    areasC.append(np.mean(aC))
    areasCplussd.append(np.mean(aC) + np.std(aC))
    areasCminussd.append(np.mean(aC) - np.std(aC))
    sizesC.append(np.mean(sC))
    curvC.append(np.mean(cC))
    curvCplussd.append(np.mean(cC) + np.std(cC))
    curvCminussd.append(np.mean(cC) - np.std(cC))
    curvCm.append(np.mean(cCm))
    curvCmplussd.append(np.mean(cCm) + np.std(cCm))
    curvCmminussd.append(np.mean(cCm) - np.std(cCm))
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

fig,ax = plt.subplots(figsize = (12, 10), dpi=200)
#plt.title('Mean areas')
ax.plot(realtime, areasC, c="#4f81c2", label="Control", linewidth=6.0)
ax.fill_between(realtime, areasCplussd, areasCminussd, color='#4f81c2', alpha=.2)
ax.plot(realtime,areasI, c="#dc756d", label="Stimulated", linewidth=6.0)
ax.fill_between(realtime, areasIplussd, areasIminussd, color='#dc756d', alpha=.2)
ax.set_xlabel('Time (h)')
ax.set_ylabel('Relative projected colony area')
ax.legend(loc = "lower left", fontsize = 40)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show()
fig.savefig(FolderToSave + "/" + d1+  'MeanAreasfig.png', bbox_inches='tight')


#Okay add graphs for main and minor axis separately. Then do of all axis together? 
fig,ax = plt.subplots(figsize = (12, 10), dpi=200)
#plt.title('Mean curvatures Main axis')
ax.plot(realtime, curvC, c="#4f81c2", label="Control", linewidth=6.0)
ax.fill_between(realtime, curvCplussd, curvCminussd, color='#4f81c2', alpha=.2)
ax.plot(realtime,curvI, c="#dc756d", label="Stimulated", linewidth=6.0)
ax.fill_between(realtime, curvIplussd, curvIminussd, color='#dc756d', alpha=.2)
ax.set_xlabel('Time (h)')
ax.set_ylabel('Average curvature')
ax.legend(loc = "upper left", fontsize = 40)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show()
plt.savefig(FolderToSave + "/" + d1 + 'MeanCurvaturesMainfig.png', bbox_inches='tight')

# Now for supplementary figure, I will make a graph with all the original curvatures considered
plt.rcParams['xtick.labelsize']=45
plt.rcParams['ytick.labelsize']=25
#Plotting controls
plt.figure(figsize = (90, 10), dpi=200)
for i in range(len(listofcurvC)):
    ax = plt.subplot(1,7,i+1)
    plt.plot(realtime, listofcurvC[i], c="#4f81c2", linewidth=8.0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('Time (h)')
    plt.ylabel('Curvature')
    plt.ylim([-0.015,0.02])
#plt.show()
plt.savefig(FolderToSave + "/" + d1 + 'AllCurvaturesControl.png', bbox_inches='tight')

#Ploting stimulated
plt.figure(figsize = (90, 10), dpi=200)
for i in range(len(listofcurvC)):
    ax = plt.subplot(1,7,i+1)
    plt.plot(realtime, listofcurvI[i], c="#dc756d", linewidth=8.0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('Time (h)')
    plt.ylabel('Curvature')
    plt.ylim([-0.015,0.02])
#plt.show()
plt.savefig(FolderToSave + "/" + d1 + 'AllCurvaturesStim.png', bbox_inches='tight')

#We need to export this data for source data file. 
import pandas as pd
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"

dataexp = pd.DataFrame(list(zip(realtime, listofcurvC[0], listofcurvC[1], listofcurvC[2], listofcurvC[3], listofcurvC[4], listofcurvC[5], listofcurvC[6])), 
                       columns = ["Time (h)", "curvature colony 1", "curvature colony 2", "curvature colony 3", "curvature colony 4", "curvature colony 5", "curvature colony 6", "curvature colony 7"])

dataexp.to_csv(filepath + '\\controlcoloniescurvatureall.csv', index = False)

dataexp = pd.DataFrame(list(zip(realtime, listofcurvI[0], listofcurvI[1], listofcurvI[2], listofcurvI[3], listofcurvI[4], listofcurvI[5], listofcurvI[6])), 
                       columns = ["Time (h)", "curvature colony 1", "curvature colony 2", "curvature colony 3", "curvature colony 4", "curvature colony 5", "curvature colony 6", "curvature colony 7"])

dataexp.to_csv(filepath + '\\stimulatedcoloniescurvatureall.csv', index = False)

#C450V controls:
DataOutputFolderC450V = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\OptoShroom Folding\C450V data"
    

listofnamesC450V = listdir(DataOutputFolderC450V)
print(listofnamesC450V)
NumbColoniesC450V = len(listofnamesC450V)
print(NumbColoniesC450V)

listoffilesC450V = []
for i in range(len(listofnamesC450V)):
    loadfileC450V = join(DataOutputFolderC450V, listofnamesC450V[i])
    temp = np.load(loadfileC450V, allow_pickle=True)
    listoffilesC450V.append(temp)

# Cool! from here we will have all the data to play with! 
print(listoffilesC450V[0][0])

time = 13 
#Area, ax size and curvature averages of main axis, clasiffied by C or Ind
listofareasCC450V = []
listofareasIC450V = []
listofsizesCC450V = []
listofsizesIC450V = []
listofcurvCC450V = []
listofcurvIC450V = []
listofcurvCmC450V = [] # for minor axis
listofcurvImC450V = [] #for minor axis
for i in range(NumbColoniesC450V):
    if listoffiles[i][0] == "C":
        listofareasCC450V.append(listoffilesC450V[i][5])
        listofsizesCC450V.append(listoffilesC450V[i][2])
        listofcurvCC450V.append(listoffilesC450V[i][6])
        listofcurvCmC450V.append(listoffilesC450V[i][7])
    else:
        listofareasIC450V.append(listoffilesC450V[i][5])
        listofsizesIC450V.append(listoffilesC450V[i][2])
        listofcurvIC450V.append(listoffilesC450V[i][6])
        listofcurvImC450V.append(listoffilesC450V[i][7])

#Cool now let's do the averages and sd:
areasCC450V = []
areasCplussdC450V = []
areasCminussdC450V = []
areasIC450V = []
areasIplussdC450V = []
areasIminussdC450V = []
sizesCC450V = []
sizesIC450V = []
curvCC450V = []
curvCplussdC450V = []
curvCminussdC450V = []
curvCmC450V = []
curvCmplussdC450V = []
curvCmminussdC450V = []
curvIC450V = []
curvIplussdC450V = []
curvIminussdC450V = []
curvImC450V = []
curvImplussdC450V = []
curvImminussdC450V = []
for t in range(time):
    aC =[]
    aI =[]
    sC =[]
    sI =[]
    cC =[]
    cI =[]
    cCm = []
    cIm = []
    for i in range(len(listofareasCC450V)):
        aC.append(listofareasCC450V[i][t]/listofareasCC450V[i][0]) #Normalizing data!
        sC.append(listofsizesCC450V[i][t]/listofsizesCC450V[i][0]) #Normalizing data!
        cC.append(listofcurvCC450V[i][t]) #Not for curvature!
        cCm.append(listofcurvCmC450V[i][t])
        aI.append(listofareasIC450V[i][t]/listofareasIC450V[i][0]) #Normalizing data!
        sI.append(listofsizesIC450V[i][t]/listofsizesIC450V[i][0]) #Normalizing data!
        cI.append(listofcurvIC450V[i][t]) #Not for curvature!
        cIm.append(listofcurvImC450V[i][t]) 
    areasCC450V.append(np.mean(aC))
    areasCplussdC450V.append(np.mean(aC) + np.std(aC))
    areasCminussdC450V.append(np.mean(aC) - np.std(aC))
    sizesCC450V.append(np.mean(sC))
    curvCC450V.append(np.mean(cC))
    curvCplussdC450V.append(np.mean(cC) + np.std(cC))
    curvCminussdC450V.append(np.mean(cC) - np.std(cC))
    curvCmC450V.append(np.mean(cCm))
    curvCmplussdC450V.append(np.mean(cCm) + np.std(cCm))
    curvCmminussdC450V.append(np.mean(cCm) - np.std(cCm))
    areasIC450V.append(np.mean(aI))
    areasIplussdC450V.append(np.mean(aI) + np.std(aI))
    areasIminussdC450V.append(np.mean(aI) - np.std(aI))
    sizesIC450V.append(np.mean(sI))
    curvIC450V.append(np.mean(cI))
    curvIplussdC450V.append(np.mean(cI) + np.std(cI))
    curvIminussdC450V.append(np.mean(cI) - np.std(cI))
    curvImC450V.append(np.mean(cIm))
    curvImplussdC450V.append(np.mean(cIm) + np.std(cIm))
    curvImminussdC450V.append(np.mean(cIm) - np.std(cIm))

#Plot all curvatures for supp figure> 
#Ploting stimulated
plt.figure(figsize = (90, 10), dpi=200)
for i in range(len(listofcurvCC450V)):
    ax = plt.subplot(1,7,i+1)
    plt.plot(realtime, listofcurvIC450V[i], c="#60D394", linewidth=8.0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('Time (h)')
    plt.ylabel('Curvature')
    plt.ylim([-0.015,0.02])
#plt.show()
plt.savefig(FolderToSave + "/" + d1 + 'AllCurvaturesStimC450V.png', bbox_inches='tight')

#Also export this data:
dataexp = pd.DataFrame(list(zip(realtime, listofcurvIC450V[0], listofcurvIC450V[1], listofcurvIC450V[2], listofcurvIC450V[3], listofcurvIC450V[4])), 
                       columns = ["Time (h)", "curvature colony 1", "curvature colony 2", "curvature colony 3", "curvature colony 4", "curvature colony 5",])

dataexp.to_csv(filepath + '\\C450Vcoloniescurvatureall.csv', index = False)


plt.rcParams.update({'font.size': 45})
plt.rcParams['ytick.labelsize']=45
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
#Now we can add the replot:
fig,ax = plt.subplots(figsize = (14, 14), dpi=200)
#plt.title('Mean areas')
ax.plot(realtime,areasIC450V, c="#60D394", label="Stimulated C450V mutant", linewidth=6.0)
ax.fill_between(realtime, areasIplussdC450V, areasIminussdC450V, color='#60D394', alpha=.2)
ax.plot(realtime, areasC, c="#4f81c2", label="Non-stimulated", linewidth=6.0)
ax.fill_between(realtime, areasCplussd, areasCminussd, color='#4f81c2', alpha=.2)
ax.plot(realtime,areasI, c="#dc756d", label="Stimulated", linewidth=6.0)
ax.fill_between(realtime, areasIplussd, areasIminussd, color='#dc756d', alpha=.2)
ax.set_xlabel('Time (h)')
ax.set_ylabel('Relative projected colony area')
ax.legend(loc = "lower left", fontsize = 40)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show()
fig.savefig(FolderToSave + "/" + 'MeanAreasfigfull.png', bbox_inches='tight')


#Data export:
dataexp = pd.DataFrame(list(zip(realtime, areasIC450V, areasIplussdC450V, areasIminussdC450V, areasC, areasCplussd, areasCminussd, areasI, areasIplussd, areasIminussd,)), 
                       columns = ["Time (h)", "area C450V", "area + std C450V", "area - std C450V", "area non-stimulated", "area + std non-stimulated", "area - std non-stimulated", "area stimulated", "area + std stimulated", "area - std stimulated"])

dataexp.to_csv(filepath + '\\foldingareas.csv', index = False)


#Now for curvatures!
fig,ax = plt.subplots(figsize = (14, 14), dpi=200)
#plt.title('Mean curvatures Main axis')
ax.plot(realtime,curvIC450V, c="#60D394", label="Stimulated C450V mutant", linewidth=6.0)
ax.fill_between(realtime, curvIplussdC450V, curvIminussdC450V, color='#60D394', alpha=.2)
ax.plot(realtime, curvC, c="#4f81c2", label="Non-stimulated", linewidth=6.0)
ax.fill_between(realtime, curvCplussd, curvCminussd, color='#4f81c2', alpha=.2)
ax.plot(realtime,curvI, c="#dc756d", label="Stimulated", linewidth=6.0)
ax.fill_between(realtime, curvIplussd, curvIminussd, color='#dc756d', alpha=.2)
ax.set_xlabel('Time (h)')
ax.set_ylabel('Average curvature ' + r'$ (μm^{-1})$')
ax.legend(loc = "upper left", fontsize = 40)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show()
fig.savefig(FolderToSave + "/" + 'MeanCurvaturesMainfigfull.png', bbox_inches='tight')

#Data export:
dataexp = pd.DataFrame(list(zip(realtime, curvIC450V, curvIplussdC450V, curvIminussdC450V, curvC, curvCplussd, curvCminussd, curvI, curvIplussd, curvIminussd,)), 
                       columns = ["Time (h)", "curvature C450V", "curvature + std C450V", "curvature - std C450V", "curvature non-stimulated", "curvature + std non-stimulated", "curvature - std non-stimulated", "curvature stimulated", "curvature + std stimulated", "curvature - std stimulated"])

dataexp.to_csv(filepath + '\\foldingcurv.csv', index = False)
