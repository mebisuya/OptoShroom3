# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 14:00:57 2020

@author: ara
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import date
today = date.today()
d1 = today.strftime("%Y_%m_%d")

from os.path import join
Folderday5 = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Day 5\Measurements"
Folderday6 = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Day 6\Measurements"
Folderday7 = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Day 7\Measurements"
Folderday8 = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Day 8\Measurements"
Folderday11 = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\EYECUP Organoid DATA\Day 11\Measurements"
FolderToSave = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\Eyecup Constriction Analysis\Graph"
from os import listdir
from numpy import genfromtxt
###IMPORTING DATA
#DAY 5
listofnames5 = listdir(Folderday5)
print(listofnames5)
NumbSamples5 = len(listofnames5)
print(NumbSamples5)

filesday5 = np.zeros((0,6))
for i in range(len(listofnames5)):
    loadfile = join(Folderday5, listofnames5[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesday5 = np.vstack((filesday5, temparray[0]))
print(filesday5) #3 first obs are time0, 3 later are time 12

#DAY 6
listofnames6 = listdir(Folderday6)
print(listofnames6)
NumbSamples6 = len(listofnames6)
print(NumbSamples6)

from numpy import genfromtxt
filesday6 = np.zeros((0,6))
for i in range(len(listofnames6)):
    loadfile = join(Folderday6, listofnames6[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesday6 = np.vstack((filesday6, temparray[0]))
print(filesday6) #3 first obs are time0, 3 later are time 12

#DAY 7
listofnames7 = listdir(Folderday7)
print(listofnames7)
NumbSamples7 = len(listofnames7)
print(NumbSamples7)

filesday7 = np.zeros((0,6))
for i in range(len(listofnames7)):
    loadfile = join(Folderday7, listofnames7[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesday7 = np.vstack((filesday7, temparray[0]))
print(filesday7) #3 first obs are time0, 3 later are time 12

#DAY 8
listofnames8 = listdir(Folderday8)
print(listofnames8)
NumbSamples8 = len(listofnames8)
print(NumbSamples8)

filesday8 = np.zeros((0,6))
for i in range(len(listofnames8)):
    loadfile = join(Folderday8, listofnames8[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesday8 = np.vstack((filesday8, temparray[0]))
print(filesday8) #3 first obs are time0, 3 later are time 12

#DAY 11
listofnames11 = listdir(Folderday11)
print(listofnames11)
NumbSamples11 = len(listofnames11)
print(NumbSamples11)

filesday11 = np.zeros((0,6))
for i in range(len(listofnames11)):
    loadfile = join(Folderday11, listofnames11[i])
    temp = genfromtxt(loadfile, delimiter=',')
    #we only want length measurements! 
    temparray = np.zeros((1,0))
    for n in range(1,7):
        temparray = np.column_stack((temparray,temp[n,6]))
    filesday11 = np.vstack((filesday11, temparray[0]))
print(filesday11) #3 first obs are time0, 3 later are time 12

#Now on to data processing, then making graphs
Day5Sizes = []
Day5end = []
Day5Diffum = []
Day5Diffper = []
for i in range(NumbSamples5):
    tempsize = (filesday5[i,0] + filesday5[i,1] + filesday5[i,2])/3
    tempsizeend = (filesday5[i,3] + filesday5[i,4] + filesday5[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    Day5Sizes.append(tempsize)
    Day5Diffum.append(tempdiff)
    Day5Diffper.append(tempper)
    Day5end.append(tempsizeend)
print(Day5Sizes)
print(Day5Diffum)
print(Day5Diffper)

#Day6
Day6Sizes = []
Day6end = []
Day6Diffum = []
Day6Diffper = []
for i in range(NumbSamples6):
    tempsize = (filesday6[i,0] + filesday6[i,1] + filesday6[i,2])/3
    tempsizeend = (filesday6[i,3] + filesday6[i,4] + filesday6[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    Day6Sizes.append(tempsize)
    Day6Diffum.append(tempdiff)
    Day6Diffper.append(tempper)
    Day6end.append(tempsizeend)
print(Day6Sizes)
print(Day6Diffum)
print(Day6Diffper)

#Day7
Day7Sizes = []
Day7end = []
Day7Diffum = []
Day7Diffper = []
for i in range(NumbSamples7):
    tempsize = (filesday7[i,0] + filesday7[i,1] + filesday7[i,2])/3
    tempsizeend = (filesday7[i,3] + filesday7[i,4] + filesday7[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    Day7Sizes.append(tempsize)
    Day7Diffum.append(tempdiff)
    Day7Diffper.append(tempper)
    Day7end.append(tempsizeend)
print(Day7Sizes)
print(Day7Diffum)
print(Day7Diffper)

#Day8
Day8Sizes = []
Day8end = []
Day8Diffum = []
Day8Diffper = []
for i in range(NumbSamples8):
    tempsize = (filesday8[i,0] + filesday8[i,1] + filesday8[i,2])/3
    tempsizeend = (filesday8[i,3] + filesday8[i,4] + filesday8[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    Day8Sizes.append(tempsize)
    Day8Diffum.append(tempdiff)
    Day8Diffper.append(tempper)
    Day8end.append(tempsizeend)
print(Day8Sizes)
print(Day8Diffum)
print(Day8Diffper)

#Day11
Day11Sizes = []
Day11end = []
Day11Diffum = []
Day11Diffper = []
for i in range(NumbSamples11):
    tempsize = (filesday11[i,0] + filesday11[i,1] + filesday11[i,2])/3
    tempsizeend = (filesday11[i,3] + filesday11[i,4] + filesday11[i,5])/3
    tempdiff = tempsizeend- tempsize
    tempper = tempdiff/tempsize #Percentage of increase
    Day11Sizes.append(tempsize)
    Day11Diffum.append(tempdiff)
    Day11Diffper.append(tempper)
    Day11end.append(tempsizeend)
print(Day11Sizes)
print(Day11Diffum)
print(Day11Diffper)

##Go make those graphs
initialsizes = [Day5Sizes, Day6Sizes, Day7Sizes, Day8Sizes, Day11Sizes]
endsizes = [Day5end, Day6end, Day7end, Day8end, Day11end]
sizediff = [Day5Diffum, Day6Diffum, Day7Diffum, Day8Diffum, Day11Diffum]
perdiff = [Day5Diffper, Day6Diffper, Day7Diffper, Day8Diffper, Day11Diffper]

#Now running paired t-test with sampels from day 5-8
import math

initialsizest = Day5Sizes + Day6Sizes + Day7Sizes + Day8Sizes
endsizest = Day5end + Day6end + Day7end + Day8end
print(endsizest)
nobs = len(initialsizest) #Number of observations

#Cool let's check for normal distribution and then do it with scipy stats
#Will do saphiro-wilk test for normality although histograms should work too

from scipy.stats import shapiro, ttest_rel
# seed the random number generator

# normality test
stat, p = shapiro(initialsizest)
print('Statistics=%.3f, p=%.3f' % (stat, p))
# interpret
alpha = 0.05
if p > alpha:
	print('Sample looks Gaussian (fail to reject H0)')
else:
	print('Sample does not look Gaussian (reject H0)')
    
stat, p = shapiro(endsizest)
print('Statistics=%.3f, p=%.3f' % (stat, p))
# interpret
alpha = 0.05
if p > alpha:
	print('Sample looks Gaussian (fail to reject H0)')
else:
	print('Sample does not look Gaussian (reject H0)')

stat, p = shapiro(sizedifft)
print('Statistics=%.3f, p=%.3f' % (stat, p))
# interpret
alpha = 0.05
if p > alpha:
	print('Sample looks Gaussian (fail to reject H0)')
else:
	print('Sample does not look Gaussian (reject H0)')
    
##All looks reasonably gaussian!
#Now finally onto running the test: 
tstat, pval = ttest_rel(initialsizest, endsizest)
#Ok P value is very small! 0.0000000004821
print(pval)

#I need to simplify data: 

plt.rcParams.update({'font.size': 16})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


fig, ax = plt.subplots(figsize = (2, 6), dpi=200)
#plt.title('Size variation in stimulated samples')
ax.scatter([1 for i in initialsizest], initialsizest, c = '#4f81c2', alpha =0.4)
ax.scatter([2 for i in endsizest], endsizest, c = '#dc756d', alpha =0.4)
for i in range(len(initialsizest)):
    ax.plot([1,2],[initialsizest[i], endsizest[i]], c='black', alpha =0.4)
ax.set_ylabel('NE thickness (' + r'$\mu$' + 'm)')
ax.set_xticks([1,2], ['0 min', '55 min'])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.text(1, max(endsizest)+5, "p = " + "%.2e"%pval )
ax.plot([1,2],[max(endsizest)+3,max(endsizest)+3], c='k')

fig.savefig(FolderToSave + '/' + 'EyecupThicknessChange.png', bbox_inches='tight')
plt.show()


#Exporting source data:
import pandas as pd
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"
dataexp = pd.DataFrame(list(zip(initialsizest, endsizest)), columns = ["Thickness 0 min", "Thickness 55 min"])

dataexp.to_csv(filepath + '\\eyecupthickness.csv', index = False)
