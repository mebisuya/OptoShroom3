# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 08:56:20 2021

@author: ara
"""

#PIV averager version 2
#This time we will first select vectors by area and average them, then find the 
#average between experiments, this way the SD will show differences between exp. 

#First we need to import all data! 
import numpy as np
import math
import matplotlib.pyplot as plt
from os import listdir
from os.path import join
from numpy import genfromtxt, savetxt

#Make plot font bigger!!
plt.rcParams.update({'font.size': 22})

#The point is that for every experiment we have a txt file for every timepoint with
#all vector values for that timepoint.
#Need to: Import all timepoints for each experiment

Folderpath = r"C:/Users/ara/Documents/Python Scripts/PIV analysis constrictions/PIV Output/Single Square stimulations/Samples"
FolderToSaveGraphs = r"C:/Users/ara/Documents/Python Scripts/PIV analysis constrictions/PIV Output/Single Square stimulations/Outputs"
graphnames = 'NewAverage' 
FolderToSaveData = r"C:\Users\ara\Documents\Python Scripts\PIV analysis constrictions\PIV Output\Single Square stimulations\Averaged data"
pixelsize = 0.276 #In um Should always be the same! 

listofnames = listdir(Folderpath)
print(listofnames)
snumb = len(listofnames) #Numberof samples
print(snumb)

files = []
for i in range(len(listofnames)):
    templistfiles = listdir(join(Folderpath, listofnames[i]))
    tempfile = []
    for b in range(len(templistfiles)):
        loadfile = join(Folderpath, listofnames[i], templistfiles[b])
        temp = genfromtxt(loadfile, delimiter=' ')
        tempfile.append(temp[:,0:5])#It is read as a NP array!
    files.append(tempfile)
print(files[0][0]) # Order is experiment, timepoint, vector. 
tps = len(files[3]) #Timepoints
print(tps)
nvec = len(files[0][0])
print(nvec)

#Objective, to make a function that takes the file[0] numpy array which is [196,5] [xpos,ypos, xcomp, ycomp, magn]
#And get from it another numpy array with xpos, ypos, radialcomponentwithsign] in which we will apply filters to plot. 

#First we need vector projection function: 
#Conceptually with v=(X,Y) u=(W,Z) Projection of v onto u.

#x value ((x*w + y*z)/(abs(w+z)**2))*w 
#y value ((x*w + y*z)/(abs(w+z)**2))*z

#In a function: 

def vproj(x,y,w,z): #projection of first over second
  projx = ((x * w + y * z)/((w**2 + z**2))) * w
  projy = ((x*w + y*z)/((w**2+z**2)))*z
  vect = [projx, projy]
  return vect


def projector(array):
    endarray = np.zeros((0,5))
    for i in range(array.shape[0]):
        #For each vector
        posvect = [256 - array[i,0], 256 -array[i,1]] #positional vector to reference point (center)
        pivvect = [array[i,2], array[i,3]] #piv vector
        projvect = vproj(pivvect[0],pivvect[1], posvect[0], posvect[1])
        #Quick way I thought to figure out sign:
        posmag = math.sqrt(posvect[0]**2 + posvect[1]**2) #  magnitude of positional vector == distance to center!
        sumvect = [posvect[0] + projvect[0], posvect[1] + projvect[1]]
        summag = math.sqrt(sumvect[0]**2 + sumvect[1]**2) #  magnitude of positional vector
        if posmag > summag:#Means that positional vector and projection have same direction! 
            mag = math.sqrt(projvect[0]**2 + projvect[1]**2) *-1# it is the magnitude in the direction to the center.
        else: #if equal or smaller
            mag = math.sqrt(projvect[0]**2 + projvect[1]**2)  # It goes in the opposite direction!
        finalset = [array[i,0], array[i,1], mag, array[i,4], posmag] #Let's include initial magnitude to compare with
        endarray = np.vstack((endarray, finalset))
    return endarray #This end array is, [xposition, yposition, radialmagnitude, total magnitude, distance to origin]

test = projector(files[0][2])
test2 = projector(files[0][7])
test3 = projector(files[0][16])

print("test 1 total mean is " + str(np.mean(test[:,3])) + " projected mean is " + str(np.mean(test[:,2])))
print("test 2 total mean is " + str(np.mean(test2[:,3])) + " projected mean is " + str(np.mean(test2[:,2])))
print("test 3 total mean is " + str(np.mean(test3[:,3])) + " projected mean is " + str(np.mean(test3[:,2])))

#Cool that seems to work! We will do it for all files
#Then we need to filter the vectors to sort them into defined areas

def selector(endarray):
    outerarray = np.zeros((0,5))
    middlearray = np.zeros((0,5))
    innerarray = np.zeros((0,5))
    for i in range(endarray.shape[0]):
        if endarray[i,0]>192 and endarray[i,1]>192 and endarray[i,0]<320 and endarray[i,1]<320:
            innerarray = np.vstack((innerarray, endarray[i,:]))
        elif endarray[i,0]>112 and endarray[i,1]>112 and endarray[i,0]<400 and endarray[i,1]<400:
            middlearray = np.vstack((middlearray, endarray[i,:]))
        else:
            outerarray = np.vstack((outerarray, endarray[i,:]))
    return outerarray, middlearray, innerarray

#Now I will do the same processing i did for average vectors but for all experiments
nexp = len(files) #Number of experiments

allvec = [] #List with all vectors and all timepoints
allmeanproj = []
outmeanproj = []
midmeanproj = []
inmeanproj = []
for i in range(nexp):
    allvect = [] #List with all vectors and all timepoints
    allmeanprojt = []
    outmeanprojt = []
    midmeanprojt = []
    inmeanprojt = []
    for n in range(tps):
        print(n)
        tempall = projector(files[i][n])
        tempout, tempmid, tempin = selector(tempall)
        
        allvect.append(tempall)
        #Let's add means
        allmeanprojt.append(np.mean(tempall[:,2])*pixelsize)
        outmeanprojt.append(np.mean(tempout[:,2])*pixelsize)
        midmeanprojt.append(np.mean(tempmid[:,2])*pixelsize)
        inmeanprojt.append(np.mean(tempin[:,2])*pixelsize)
    allvec.append(allvect)
    allmeanproj.append(allmeanprojt)
    outmeanproj.append(outmeanprojt)
    midmeanproj.append(midmeanprojt)
    inmeanproj.append(inmeanprojt)

#Ok, here we have mean velocity values for each per each timepoint in each experiment
#Now we want to average by experiment and get SDs


allvelmean = []
allvelsd = []
outervelmean = []
outervelsd = []
midvelmean = []
midvelsd = []
invelmean = []
invelsd = []
for n in range(tps):
    allvelt = []
    outervelt = []
    midvelt = []
    invelt = []
    for i in range(nexp):
        allvelt.append(allmeanproj[i][n]/2) #No division by 2 we sum up total displacements
        outervelt.append(outmeanproj[i][n]/2)
        midvelt.append(midmeanproj[i][n]/2)
        invelt.append(inmeanproj[i][n]/2)
    allvelmean.append(np.mean(allvelt))
    allvelsd.append(np.std(allvelt))
    outervelmean.append(np.mean(outervelt))
    outervelsd.append(np.std(outervelt))
    midvelmean.append(np.mean(midvelt))
    midvelsd.append(np.std(midvelt))
    invelmean.append(np.mean(invelt))
    invelsd.append(np.std(invelt))

print(midvelmean[5])
print(midvelsd[5])

#Results are differentwhen changing order of averaging?
plt.rcParams.update({'font.size': 7})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['pdf.fonttype'] = 42

fig, ax = plt.subplots(figsize = (2.5, 2))
#plt.title('Classified measured velocities (radial')
for i in range(8, 68):
    ax.axvspan(i + 0.333, i + 0.917, color = '#a4c3f6', alpha = 0.3, fill = "True", lw= 0)
ax.plot([2*i for i in range(tps)], outervelmean, linewidth = 1, c='g', label= 'Outer area', alpha =0.8)
ax.plot([2*i for i in range(tps)], midvelmean, linewidth = 1, color='#4f81c2', label= 'Adjacent area', alpha =0.8)
ax.plot([2*i for i in range(tps)], invelmean, linewidth = 1, c='#dc756d', label= 'Stimulated area', alpha =0.8)
ax.set_ylabel('Radial velocity (μm/min)')
ax.set_xlabel('Time (min)')
ax.legend(fontsize = 6)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.show()
plt.savefig(FolderToSaveGraphs + '/' + graphnames + 'ClassVectorsSpeeds.pdf', bbox_inches='tight', transparent=True)

#Great! Let's export the source data:
import pandas as pd
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"

dataexp = pd.DataFrame(list(zip([2*i for i in range(tps)], outervelmean, midvelmean, invelmean)),
                       columns = ["Time (min)", "Outer area", "Adjacent area", "Stimulated area"])

dataexp.to_csv(filepath + '\\PIVaverages.csv', index = False)

#Now for accumulated displacement! 



#first, make new list summing

sumallmeanproj = []
sumoutmeanproj = []
summidmeanproj = []
suminmeanproj = []
for i in range(nexp):
    allmeanprojt = []
    outmeanprojt = []
    midmeanprojt = []
    inmeanprojt = []
    for n in range(tps):
        if n == 0:
            allmeanprojt.append((allmeanproj[i][n]))
            outmeanprojt.append((outmeanproj[i][n]))
            midmeanprojt.append((midmeanproj[i][n]))
            inmeanprojt.append((inmeanproj[i][n]))
        else:
            allmeanprojt.append((allmeanprojt[n-1] + allmeanproj[i][n]))
            outmeanprojt.append((outmeanprojt[n-1] + outmeanproj[i][n]))
            midmeanprojt.append((midmeanprojt[n-1] + midmeanproj[i][n]))
            inmeanprojt.append((inmeanprojt[n-1] + inmeanproj[i][n]))
    sumallmeanproj.append(allmeanprojt)
    sumoutmeanproj.append(outmeanprojt)
    summidmeanproj.append(midmeanprojt)
    suminmeanproj.append(inmeanprojt)                    

print(summidmeanproj[2])

#Now average experiments

allsummean = []
allsumsd = []
outersummean = []
outersumsd = []
midsummean = []
midsumsd = []
insummean = []
insumsd = []
for n in range(tps):
    allsumt = []
    outersumt = []
    midsumt = []
    insumt = []
    for i in range(nexp):
        allsumt.append(sumallmeanproj[i][n]) #2 is the number of minutes between frames
        outersumt.append(sumoutmeanproj[i][n])
        midsumt.append(summidmeanproj[i][n])
        insumt.append(suminmeanproj[i][n])
    allsummean.append(np.mean(allsumt))
    allsumsd.append(np.std(allsumt))
    outersummean.append(np.mean(outersumt))
    outersumsd.append(np.std(outersumt))
    midsummean.append(np.mean(midsumt))
    midsumsd.append(np.std(midsumt))
    insummean.append(np.mean(insumt))
    insumsd.append(np.std(insumt))


fig, ax = plt.subplots(figsize = (2.5, 2))
#plt.title('Classidied radial displacements')
for i in range(8, 68):
    ax.axvspan(i + 0.333, i + 0.917, color = '#a4c3f6', alpha = 0.3, fill = "True", lw= 0)
ax.plot([2*i for i in range(tps)], outersummean, linewidth = 1, c='g', label= 'Outer area', alpha = 0.8)
ax.plot([2*i for i in range(tps)], midsummean, linewidth = 1, color='#4f81c2', label= 'Adjacent area', alpha = 0.8)
ax.plot([2*i for i in range(tps)], insummean, linewidth = 1, color='#dc756d', label= 'Stimulated area', alpha = 0.8)
ax.set_ylabel('Radial displacement (μm)')
ax.set_xlabel('Time (min)')
ax.legend(fontsize = 6, loc = "lower center")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.show()
plt.savefig(FolderToSaveGraphs + '/' + graphnames + 'ClassVectorsDisp.pdf', bbox_inches='tight', transparent=True)


#Exporting data
dataexp = pd.DataFrame(list(zip([2*i for i in range(tps)], outersummean, midsummean, insummean)),
                       columns = ["Time (min)", "Outer area", "Adjacent area", "Stimulated area"])

dataexp.to_csv(filepath + '\\PIVsums.csv', index = False)


print(midsummean[-5])
print(midsumsd[-5])
#Outer area
print(outersummean[-5])
print(outersumsd[-5])
#Inner area
print(insummean[-5])
print(insumsd[-5])

plt.figure(figsize = (12, 8), dpi=200)
#plt.title('Classified measured velocities (radial')
for i in range(8, 68):
    plt.axvspan(i + 0.333, i + 0.917, color = '#a4c3f6', alpha = 0.3)
#plt.plot([2*i for i in range(tps)], [outervelmean[n] for n in range(tps)], linewidth = 5, c='g', label= 'Outer region', alpha =0.8)
plt.plot([2*i for i in range(tps)], [summidmeanproj[0][n] for n in range (tps)], linewidth = 5, color='#4f81c2', label= 'Adjacent region experiments', alpha =0.8)
plt.plot([2*i for i in range(tps)], [summidmeanproj[1][n] for n in range (tps)], linewidth = 5, color='#4f81c2', alpha =0.8)
plt.plot([2*i for i in range(tps)], [summidmeanproj[2][n] for n in range (tps)], linewidth = 5, color='#4f81c2', alpha =0.8)
plt.plot([2*i for i in range(tps)], [summidmeanproj[3][n] for n in range (tps)], linewidth = 5, color='#4f81c2', alpha =0.8)
plt.plot([2*i for i in range(tps)], [summidmeanproj[4][n] for n in range (tps)], linewidth = 5, color='#4f81c2', alpha =0.8)
plt.plot([2*i for i in range(tps)], [summidmeanproj[5][n] for n in range (tps)], linewidth = 5, color='#4f81c2', alpha =0.8)
plt.plot([2*i for i in range(tps)], [midsummean[n] for n in range (tps)], linewidth = 5, color="#dc756d", label= 'Adjacent region mean', alpha =0.8)
#plt.plot([2*i for i in range(tps)], [invelmean[n] for n in range(tps)], linewidth = 5, c='#dc756d', label= 'Stimulated region', alpha =0.8)
plt.ylabel('Radial displacement (μm)')
plt.xlabel('Time (min)')
plt.legend(fontsize = 22)
ax = plt.axes()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.show()
plt.savefig(FolderToSaveGraphs + '/' + graphnames + 'CompositeSummStimulations.png', bbox_inches='tight')


#Supp figure also with SD?

linelist = []
linesdplus = []
linesdminus = []
for tp in range(tps):
    lineav = []
    sdplus = []
    sdminus = []
    distav = np.unique(allvec[0][tp][:,4])*pixelsize
    for i in np.unique(allvec[0][tp][:,4])*pixelsize:
        templist = []
        for e in range(nexp):
            for v in range(nvec):
                currentdist = allvec[e][tp][v,4] *pixelsize
                if currentdist == i:
                    templist.append(allvec[e][tp][v,2]*pixelsize/2)
        sd = np.std(templist)
        lineav.append(np.mean(templist))
        sdplus.append(np.mean(templist) + sd)
        sdminus.append(np.mean(templist) - sd)
    linelist.append(lineav)
    linesdplus.append(sdplus)
    linesdminus.append(sdminus)


plt.figure(figsize = (7, 7))
#plt.suptitle('Distribution of radial component in all vectors')
ax1 = plt.subplot(2,2,1)
plt.axvspan(0, 24.98, color = 'gray', alpha = 0.3, fill = "True", lw= 0)
#plt.title('before stimulation')
plt.plot(distav, linelist[4], c='k', label = "2 min before stimulation")
plt.fill_between(distav, linesdplus[4], linesdminus[4], color='gray', alpha=.3)
plt.ylabel('Radial velocity (' + r'$\mu$' + 'm/min)')
plt.legend(loc = "lower center", fontsize = 7)
#plt.xlabel('Distance to centre(um)')
plt.ylim(-0.27, 0.27)
plt.xlim(0, 90)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2 = plt.subplot(2,2,2)
plt.axvspan(0, 24.98, color = '#78e3d8', alpha = 0.3, fill = "True", lw= 0)
#plt.title('2 min stimulation')
plt.plot(distav, linelist[5], c='k', label = "2 min stimulation")
plt.fill_between(distav, linesdplus[5], linesdminus[5], color='gray', alpha=.3)
plt.legend(loc = "lower center", fontsize = 7)
#plt.xlabel('Distance to centre(um)')
plt.ylim(-0.27, 0.27)
plt.xlim(0, 90)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3 = plt.subplot(2,2,3)
plt.axvspan(0, 24.98, color = '#78e3d8', alpha = 0.3, fill = "True", lw= 0)
#plt.title('20 min stimulation')
plt.plot(distav, linelist[15], c='k', label = "20 min stimulation")
plt.fill_between(distav, linesdplus[15], linesdminus[15], color='gray', alpha=.3)
#plt.ylabel('velocity (um/min)')
plt.legend(loc = "lower center", fontsize = 7)
plt.ylabel('Radial velocity (' + r'$\mu$' + 'm/min)')
plt.xlabel('Distance to center (' + r'$\mu$' + 'm)')
plt.ylim(-0.27, 0.27)
plt.xlim(0, 90)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax4 = plt.subplot(2,2,4)
plt.axvspan(0, 24.98, color = 'gray', alpha = 0.3, fill = "True", lw= 0)
#plt.title('2 min after stimulation')
plt.plot(distav, linelist[35], c='k', label = "2 min after stimulation")
plt.fill_between(distav, linesdplus[35], linesdminus[35], color='gray', alpha=.3)
plt.xlabel('Distance to center (' + r'$\mu$' + 'm)')
plt.ylim(-0.27, 0.27)
plt.xlim(0, 90)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
plt.legend(loc = "upper center", fontsize = 7)
#plt.show()
plt.savefig(FolderToSaveGraphs + '/' + graphnames + 'SpatialDistribselected.pdf', bbox_inches='tight', transparent=True)

#Exporting data
dataexp = pd.DataFrame(list(zip(distav, linelist[4], linesdplus[4], linesdminus[4], linelist[5], linesdplus[5], linesdminus[5], linelist[15], linesdplus[15], linesdminus[15], linelist[35],  linesdplus[35], linesdminus[35])),
                       columns = ["Distance to center (um)", "2 min before stimulation", "2 min before stimulation + stand dev", "2 min before stimulation - stand dev", "2 min stimulation", "2 min stimulation + stand dev", "2 min stimulation - stand dev", "20 min stimulation",  "20 min stimulation + stand dev", "20 min stimulation - stand dev", "2 min after stimulation",  "2 min after stimulation + stand dev", "2 min after stimulation - stand dev",])

dataexp.to_csv(filepath + '\\PIVdistav.csv', index = False)
