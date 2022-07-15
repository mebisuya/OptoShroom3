# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 11:35:03 2021

@author: ara
"""

### graph maker for EB constriction analysis. 


import numpy as np
import matplotlib.pyplot as plt
from datetime import date
today = date.today()
d1 = today.strftime("%Y_%m_%d")
import math
from os.path import join
DataOutputFolder = r"C:\Users\ara\Documents\Python Scripts\Embryonic body constriction\Output Data"
CurvaturesFolder = r'C:\Users\ara\Documents\Python Scripts\Embryonic body constriction\curvatures only'
FolderToSave = r"C:\Users\ara\Documents\Python Scripts\Embryonic body constriction\Output graphs"
from os import listdir

#Make plot font bigger!!
#plt.rcParams.update({'font.size': 45})

listofnames = listdir(DataOutputFolder)
print(listofnames)
NumbColonies = len(listofnames)
print(NumbColonies)

listoffiles = []
for i in range(len(listofnames)):
    loadfile = join(DataOutputFolder, listofnames[i])
    temp = np.load(loadfile, allow_pickle=True)
    listoffiles.append(temp)

#Curvature
curv = []
curvplussd = []
curvminussd = []
for t in range(time):
    temp = []
    for n in range(NumbColonies):
        temp.append(listoffiles[n][4][t]) 
    mean = np.mean(temp)
    sd = np.std(temp)
    if sd>0.001:
        temp.remove(min(temp))
        print(sd)
        mean = np.mean(temp)
        sd = np.std(temp)
    else:
        pass
    curv.append(mean)
    curvplussd.append(mean+sd)
    curvminussd.append(mean-sd)


#Plot:
plt.figure(figsize = (8,8), dpi = 200)
plt.title('Curvature')
plt.axvspan(15, 70, color = '#78e3d8', alpha = 0.3)
plt.plot(realtime, curv, c='k', linewidth = 8)
plt.fill_between(realtime, curvplussd, curvminussd, color='k', alpha=.1)
plt.ylabel('Curvature')
plt.xlabel('time (min)')
#plt.xlabel('Distance to centre(um)')
#ax = plt.axes()
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
plt.show()
plt.savefig(FolderToSave + "/" + d1 + 'Curvature.png', bbox_inches='tight')


#Overall curvature does not change, which I think is expectable in a closed surface? (compensating each other out?)

#So let's go with the local curvatures:
listofnames = listdir(CurvaturesFolder)
print(listofnames)
NumbColonies = len(listofnames)
print(NumbColonies)

listoffiles = []
for i in range(len(listofnames)):
    loadfile = join(CurvaturesFolder, listofnames[i])
    temp = np.load(loadfile, allow_pickle=True)
    listoffiles.append(temp)

print(listoffiles[0][0,:,:]) #These are 3d arrays (time, measurements, degreesandcurvature)

def getAverageArray(listofarrays):
    #For 3D arrays in a nested list(different observations)
    averagearray = np.zeros_like(listofarrays[-1])
    standevup = np.zeros_like(listofarrays[-1])
    standevlow = np.zeros_like(listofarrays[-1])
    time, observations, dims = listofarrays[-1].shape # We assume they all have the same shape!
    for t in range(time):
        for o in range(observations):
            obs = []
            for i in range(len(listofarrays)):
                obs.append(listofarrays[i][t,o,1])
            averagearray[t,o,1] = np.mean(obs)
            averagearray[t,o,0] = listofarrays[i][t,o,0]
            standevup[t,o,1] = np.mean(obs) + np.std(obs)
            standevup[t,o,0] = listofarrays[i][t,o,0]
            standevlow[t,o,1] = np.mean(obs) - np.std(obs)
            standevlow[t,o,0] = listofarrays[i][t,o,0]
    return averagearray, standevup, standevlow

ave, sdup, sdlo = getAverageArray(listoffiles)

#What about changes in curvature?

#Let's first try to visualize a timeline.
#I will arbitrarily divide in ranges of 30degrees, although I could measure the average range the stimulation covers
time = ave.shape[0]
stimarea = [np.mean(ave[t,165:195,1]) for t in range(time)] #This is the range -15 - 15
adareaminus = [np.mean(ave[t,135:165,1]) for t in range(time)] # range - 25 to -55
adareaplus = [np.mean(ave[t,195:235,1]) for t in range(time)] # range 25 to 55
outerplus = [np.mean(ave[t,280:310,1]) for t in range(time)] # range 135 to 105
outerminus = [np.mean(ave[t,40:70,1]) for t in range(time)] # range -135 to -105

timeline = [i*5 for i in range(time)]
plt.figure()
plt.plot(timeline, stimarea, color = "red", label = "stimulated")
plt.plot(timeline, adareaminus, color = "green", label = "adjacent 1")
plt.plot(timeline, adareaplus, color = "green", label = "adjacent 1")
plt.plot(timeline, outerplus, color = "black", label = "outer")
plt.plot(timeline, outerminus, color = "black", label = "outer")
plt.show()

#Seems like we've got it!. Final graphs will be: 
    #ADD standard colors!
    
#1. Average curvature before and after stimulation
t1= 2
t2 = 14

#Obtained from separate code, average stimulation region:
lowerstimangle = -17.05
uppersimangle = 15.44

#To add expected curvature, we have initial areas:
print(initialareas)
#if they were a perfect circunference, they would be: pi*r**2
meanarea = np.mean(initialareas)
meanradius = math.sqrt(meanarea/math.pi)
meanexpectedcurvature = 1/meanradius


plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['xtick.labelsize']=45
plt.rcParams['ytick.labelsize']=35
plt.rcParams.update({'font.size': 45})

fig, ax =plt.subplots(figsize = (20,14), dpi = 200)
ax.plot([-180, 180], [meanexpectedcurvature, meanexpectedcurvature], color= 'k', label = 'Expected curvature', alpha = 0.4, linewidth=6.0)
ax.plot(ave[t1,:,0],ave[t1,:,1], color = "#4f81c2", alpha = 0.8, label = "Before stimulation", linewidth=6.0)
ax.fill_between(ave[t1,:,0],sdup[t1,:,1],sdlo[t1,:,1], color='#4f81c2', alpha=.1)
ax.plot(ave[t2,:,0],ave[t2,:,1], color = "#dc756d", alpha = 0.8, label = "1 h stimulation", linewidth=6.0)
ax.fill_between(ave[t2,:,0],sdup[t2,:,1],sdlo[t2,:,1], color='#dc756d', alpha=.1)
ax.set_ylabel("Curvature " +  r'$ (μm^{-1})$')
ax.set_xlabel("Orientation (degrees)")
ax.axvspan(lowerstimangle, uppersimangle, color = '#a4c3f6', alpha = 0.3)
ax.legend(loc = "lower left", fontsize = 40)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig(FolderToSave + "/" + d1 + 'maincurvaturesfig.png', bbox_inches='tight')
plt.show()

## Export:
import pandas as pd 
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"
dataexp = pd.DataFrame(list(zip(list(ave[t1,:,0]), list(ave[t1,:,1]), list(sdup[t1,:,1]), list(sdlo[t1,:,1]), list(ave[t2,:,1]), list(sdup[t2,:,1]), list(sdlo[t2,:,1]), )), 
                       columns = ["Degrees", "Curvature before stimulation", "Curvature before stim + standard deviation", "Curvature before stim - standard deviation", "Curvature 1 h stimulation", "Curvature 1 h stim + standard deviation", "Curvature 1 h stim - standard deviation"])

dataexp.to_csv(filepath + '\\neuroectcurvature.csv', index = False)


#I would like to make a gif of this graph!
from matplotlib.animation import FuncAnimation, PillowWriter  

#To save perimeter evolution
figure, ax =plt.subplots(figsize = (20,14), dpi = 200)
#figure.tight_layout(pad=0)
def update(t2):  
    ax.set_title("Curvature evolution in OptoShroom3 neuroectodermal organoids")
    ax.cla()
    if t2 > 2 and t2<15:
        plt.axvspan(lowerstimangle, uppersimangle, color = '#78e3d8', alpha = 0.3)
    ax.plot([-180, 180], [meanexpectedcurvature, meanexpectedcurvature], color= 'k', label = 'expected curvature', alpha = 0.4,  linewidth=6.0)
    ax.plot(ave[t2,:,0],ave[t2,:,1], color = "#dc756d", alpha = 0.8, linewidth = 6, label = 'measured curvature')
    ax.set_ylabel("curvature " +  r'$ (μm^{-1})$')
    ax.set_xlabel("orientation (degrees)")
    ax.set_ylim(0.004, 0.008)
    ax.text(-160, 0.0075, str(t2*5)+ " min", size = 'medium')
    ax.legend(loc = "lower left",  fontsize = 40)
    plt.fill_between(ave[t2,:,0],sdup[t2,:,1],sdlo[t2,:,1], color="#dc756d", alpha=.1)

    
ani = FuncAnimation(figure, update, frames = time)  
writer = PillowWriter(fps=10)  
ani.save(FolderToSave + "/" + "giftrial.tif", writer=writer)


#Export
dataexp = pd.DataFrame(list(zip(timeline, ave[:,:,1], sdup[:,:,1], sdlo[:,:,1])), 
                       columns = ["Time (min)", "Curvatures", "Curvatures + standard deviation", "Curvatures - standard deviation"])

dataexp.to_csv(filepath + '\\neuroectcurvatureALL.csv', index = False)
#2. Timecourse of the curvatures and their different degrees.


stimarea = [np.mean(ave[t,165:195,1]) for t in range(time)] #This is the range -15 - 15
ar135_165 = [np.mean(ave[t,135:165,1]) for t in range(time)]
ar105_135 = [np.mean(ave[t,105:135,1]) for t in range(time)]
ar75_105 = [np.mean(ave[t,75:105,1]) for t in range(time)]
ar45_75 = [np.mean(ave[t,45:75,1]) for t in range(time)]
ar15_45 = [np.mean(ave[t,15:45,1]) for t in range(time)]
ar195_225 = [np.mean(ave[t,195:225,1]) for t in range(time)]
ar225_255 = [np.mean(ave[t,225:255,1]) for t in range(time)]
ar255_285 = [np.mean(ave[t,225:285,1]) for t in range(time)]
ar285_315 = [np.mean(ave[t,285:315,1]) for t in range(time)]
ar315_345 = [np.mean(ave[t,315:345,1]) for t in range(time)]
ar345_15 = [np.mean((ave[t,345:,1], ave[t,:15,1])) for t in range(time)]
timeline = [i*5 for i in range(time)]

plt.figure()
plt.axvspan(timeline[3], timeline[15], color = '#78e3d8', alpha = 0.3)
plt.plot(timeline, stimarea, color = "red", label = "stimulated")
plt.plot(timeline, ar135_165, color = "green", label = "adjacent")
plt.plot(timeline, ar195_225, color = "green")
plt.plot(timeline, ar105_135, color = "black", label = "outer")
plt.plot(timeline, ar75_105, color = "black")
plt.plot(timeline, ar45_75, color = "black")
plt.plot(timeline, ar15_45, color = "black")
plt.plot(timeline, ar225_255, color = "black")
plt.plot(timeline, ar255_285, color = "black")
plt.plot(timeline, ar285_315, color = "black")
plt.plot(timeline, ar315_345, color = "black")
plt.plot(timeline, ar345_15, color = "black")
plt.ylabel("curvature (1/um)")
plt.xlabel("time (min)")
plt.legend()
plt.savefig(FolderToSave + "/" + d1 + 'curvatureevolution.png', bbox_inches='tight')
plt.show()

#Getting the relative variation to state it on the text:

ave[t2,170:190,0]
nonstims = np.concatenate((ave[t2,:165,1],ave[t2,195:,1]))
meanstimval = np.mean(ave[t2,170:190,1])
meannonstimval = np.mean(nonstims)
stimsd= np.std(ave[t2,170:190,1])
stimsd/meannonstimval*100
(1 - meanstimval/meannonstimval)*100

#alternatively compared with the same area before stimulation (probably better)


meanstimval = np.mean(ave[t2,170:190,1])
meannonstimval = np.mean(ave[t1,170:190,1])
stimsd= np.std(ave[t2,170:190,1])
stimsd/meannonstimval*100
(1 - meanstimval/meannonstimval)*100
