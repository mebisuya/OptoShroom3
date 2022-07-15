# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 12:51:18 2020

@author: ara
"""

#Optic vesicle graph maker, based on measurements made with fiji. 
# For each optic vesicles 3 measurements were made before stimulation and 1h after stimulation
# along the longest axis of the lumen of the vesicle


import numpy as np
import matplotlib.pyplot as plt
from datetime import date
today = date.today()
d1 = today.strftime("%Y_%m_%d")

from os.path import join
FolderC5 = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Vesicle Induction\C5\measurements"
FolderC5cont = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Vesicle Induction\C5\control measurements"
FolderC28 = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Vesicle Induction\C28\measurements"
FolderControls = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Vesicle Induction\Control 11\measurements"
FolderToSave = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Vesicle Induction\graphs"

plt.rcParams.update({'font.size': 22})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

from os import listdir
from numpy import genfromtxt
###IMPORTING DATA
#Clon 5
listofnames5 = listdir(FolderC5)
print(listofnames5)
NumbSamples5 = len(listofnames5)
print(NumbSamples5)

filesC5 = np.zeros((0,6))
for i in range(len(listofnames5)):
    loadfile = join(FolderC5, listofnames5[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesC5 = np.vstack((filesC5, temparray[0]))
print(filesC5) #3 first obs are time0, 3 later are time 12

#Clon 5 controls
listofnames5cont = listdir(FolderC5cont)
print(listofnames5cont)
NumbSamples5cont = len(listofnames5cont)
print(NumbSamples5cont)

filesC5cont = np.zeros((0,6))
for i in range(len(listofnames5cont)):
    loadfile = join(FolderC5cont, listofnames5cont[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesC5cont = np.vstack((filesC5cont, temparray[0]))
print(filesC5cont) #3 first obs are time0, 3 later are time 12

#Clon 28
listofnames28 = listdir(FolderC28)
print(listofnames28)
NumbSamples28 = len(listofnames28)
print(NumbSamples28)

filesC28 = np.zeros((0,6))
for i in range(len(listofnames28)):
    loadfile = join(FolderC28, listofnames28[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesC28 = np.vstack((filesC28, temparray[0]))
print(filesC28) #3 first obs are time0, 3 later are time 12

#Control c11
listofnamesControls = listdir(FolderControls)
print(listofnamesControls)
NumbSamplesControls = len(listofnamesControls)
print(NumbSamplesControls)

filesControls = np.zeros((0,6))
for i in range(len(listofnamesControls)):
    loadfile = join(FolderControls, listofnamesControls[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesControls = np.vstack((filesControls, temparray[0]))
print(filesControls) #3 first obs are time0, 3 later are time 12



#Making averages and calculations
C5Sizes = []
C5end = []
C5Diffum = []
C5Diffper = []
for i in range(NumbSamples5):
    tempsize = (filesC5[i,0] + filesC5[i,1] + filesC5[i,2])/3
    tempsizeend = (filesC5[i,3] + filesC5[i,4] + filesC5[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    C5Sizes.append(tempsize)
    C5Diffum.append(tempdiff)
    C5Diffper.append(tempper)
    C5end.append(tempsizeend)
print(C5Sizes)
print(C5Diffum)
print(C5Diffper)

C5contSizes = []
C5contend = []
C5contDiffum = []
C5contDiffper = []
for i in range(NumbSamples5cont):
    tempsize = (filesC5cont[i,0] + filesC5cont[i,1] + filesC5cont[i,2])/3
    tempsizeend = (filesC5cont[i,3] + filesC5cont[i,4] + filesC5cont[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    C5contSizes.append(tempsize)
    C5contDiffum.append(tempdiff)
    C5contDiffper.append(tempper)
    C5contend.append(tempsizeend)
print(C5contSizes)
print(C5contDiffum)
print(C5contDiffper)

C28Sizes = []
C28end = []
C28Diffum = []
C28Diffper = []
for i in range(NumbSamples28):
    tempsize = (filesC28[i,0] + filesC28[i,1] + filesC28[i,2])/3
    tempsizeend = (filesC28[i,3] + filesC28[i,4] + filesC28[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    C28Sizes.append(tempsize)
    C28Diffum.append(tempdiff)
    C28Diffper.append(tempper)
    C28end.append(tempsizeend)
print(C28Sizes)
print(C28Diffum)
print(C28Diffper)

ControlsSizes = []
Controlsend = []
ControlsDiffum = []
ControlsDiffper = []
for i in range(NumbSamplesControls):
    tempsize = (filesControls[i,0] + filesControls[i,1] + filesControls[i,2])/3
    tempsizeend = (filesControls[i,3] + filesControls[i,4] + filesControls[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    ControlsSizes.append(tempsize)
    ControlsDiffum.append(tempdiff)
    ControlsDiffper.append(tempper)
    Controlsend.append(tempsizeend)
print(ControlsSizes)
print(ControlsDiffum)
print(ControlsDiffper)

print(np.mean(C5Diffper))
#Let's make some graph!

plt.figure(figsize =(8, 8))
plt.title('Relative change in lumen diameter')
plt.boxplot([ControlsDiffper, C5contDiffper, C5Diffper, C28Diffper], widths = 0.6)
plt.scatter([1 for i in ControlsDiffper], ControlsDiffper)
plt.scatter([2 for i in C5contDiffper], C5contDiffper)
plt.scatter([3 for i in C5Diffper], C5Diffper)
plt.scatter([4 for i in C28Diffper], C28Diffper)
plt.ylabel('Relative size')
plt.xticks([1,2,3,4], ['Control', 'C5 control','C5', 'C2.8'])
#plt.show()
plt.savefig(FolderToSave + '/' + d1 +'Vesiclediametesall.png', bbox_inches='tight')

#Without c28
plt.figure(figsize =(8, 8))
plt.title('Relative change in lumen diameter')
plt.boxplot([ControlsDiffper, C5Diffper], widths = 0.6)
plt.scatter([1 for i in ControlsDiffper], ControlsDiffper)
plt.scatter([2 for i in C5Diffper], C5Diffper)
plt.ylabel('Relative size')
plt.xticks([1,2], ['Control','C5'])
#plt.show()
plt.savefig(FolderToSave + '/' + d1 +'Vesiclediametes.png', bbox_inches='tight')

#What about t-test for differences?

from scipy.stats import ttest_ind, ttest_rel

tstat, pval = ttest_ind(C5contDiffum, C5Diffum)
import seaborn as sns
import pandas as pd 

data = pd.DataFrame({'Non-stimulated': pd.Series(C5contDiffum), 'Stimulated': pd.Series(C5Diffum)})
print(data)

plt.rcParams.update({'font.size': 16})
plt.rcParams['xtick.labelsize']=14
# Create an array with the colors you want to use
colors = ["#4f81c2", "#dc756d"]
# Set your custom color palette
sns.set_palette(sns.color_palette(colors))

plt.figure(figsize = (3, 6), dpi=200)
ax = sns.swarmplot(data=data, size=5)
ax = sns.boxplot(data=data, boxprops={'facecolor':'None'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.ylabel('Lumen diameter change (' + r'$\mu$' + 'm)')
plt.text(-0.1, max(C5Diffum)+1.5, "p = " + "%.2e"%pval )
plt.plot([-0.2,1.2],[max(C5Diffum)+1,max(C5Diffum)+1], c='k')
#plt.show()
plt.savefig(FolderToSave + '/' +'Vesiclediametersfig.png', bbox_inches='tight')

#Export data:
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"


dataexp = pd.DataFrame({'Non-stimulated': pd.Series(C5contDiffum), 'Stimulated': pd.Series(C5Diffum)})
dataexp.to_csv(filepath + '\\vesicledifference.csv', index = False)


#Now finally onto running the test: 
tstat, pval = ttest_rel(initialsizest, endsizest)
#Ok P value is very small! 0.009
print(pval)

plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(figsize = (2, 6), dpi=200)
#plt.title('Lumen diameter in stimulated samples')
ax.scatter([1 for i in initialsizest], initialsizest, c = '#4f81c2', alpha =0.4)
#plt.scatter([1 for i in initialsizes[4]], initialsizes[4], c = 'purple', label = 'day11')
ax.scatter([2 for i in endsizest], endsizest, c = '#dc756d', alpha =0.4)
#plt.scatter([2 for i in endsizes[4]], endsizes[4], c = 'purple', label = 'day11')
for i in range(len(initialsizest)):
    ax.plot([1,2],[initialsizest[i], endsizest[i]], c='black', alpha =0.4)
#for i in range(len(initialsizes[4])):
    #plt.plot([1,2],[initialsizes[4][i], endsizes[4][i]], c='purple', alpha =0.4)
ax.set_ylabel('Lumen diameter (' + r'$\mu$' + 'm)')
ax.set_xticks([1,2], ['0 min', '55 min'])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.text(0.9, max(initialsizest)+5, "p = " + "%.2e"%pval )
ax.plot([0.9,2.1],[max(initialsizest)+3,max(initialsizest)+3], c='k')
fig.savefig(FolderToSave + '/' + 'EyecupDiameterChangeFig.png', bbox_inches='tight')
plt.show()

#Let's export the data:

dataexp = pd.DataFrame(list(zip(initialsizest, endsizest)), columns = ["Thickness 0 min", "Thickness 55 min"])

dataexp.to_csv(filepath + '\\vesiclesize.csv', index = False)