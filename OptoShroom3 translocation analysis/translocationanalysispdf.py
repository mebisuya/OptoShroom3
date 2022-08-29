# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 15:58:14 2021

@author: ara
"""

##New code for analyzing translocation dynamics measured with boxes in imageJ. 

#Trial to make sure how data is imported. 

import numpy as np
from os.path import join
from os import listdir


#We need 4 folders!

OptoShFolder = r"C:\Users\ara\Documents\Thesis\OptoShroom3 Paper\RawData and Python Scripts\Translocation analysis\OptoShroom3\Measurements"
OptoChFolder = r"C:\Users\ara\Documents\Thesis\OptoShroom3 Paper\RawData and Python Scripts\Translocation analysis\Optocherry\Measurements"
c450vC2Folder = r"C:\Users\ara\Documents\Thesis\OptoShroom3 Paper\RawData and Python Scripts\Translocation analysis\C450V Clone2\Measurements"
c450vC5Folder = r"C:\Users\ara\Documents\Thesis\OptoShroom3 Paper\RawData and Python Scripts\Translocation analysis\C450V Clone5\Measurements"
EyecupFolder = r'C:\Users\ara\Documents\Thesis\OptoShroom3 Paper\RawData and Python Scripts\Translocation analysis\EyecupTranslocation\Measurements'
#And for graphs and measurements
outputfolder = r"C:\Users\ara\Documents\Thesis\OptoShroom3 Paper\RawData and Python Scripts\Translocation analysis"

OptoShNames = listdir(OptoShFolder)

OptoShCyt = np.zeros((0,60))
OptoShJun = np.zeros((0,60))
for i in range(len(OptoShNames)):
    source = join(OptoShFolder, OptoShNames[i])
    temp = np.genfromtxt(source, delimiter = ',', skip_header = 1)
    if "cyt" in OptoShNames[i]:
        OptoShCyt = np.vstack((OptoShCyt, temp[:,2]))
    elif "jun" in OptoShNames[i]:
        OptoShJun = np.vstack((OptoShJun, temp[:,2]))
    else:
        print("variable mislabeled")
        
#Same for the other experiments!!!

OptoChNames = listdir(OptoChFolder)

OptoChCyt = np.zeros((0,60))
OptoChJun = np.zeros((0,60))
for i in range(len(OptoChNames)):
    source = join(OptoChFolder, OptoChNames[i])
    temp = np.genfromtxt(source, delimiter = ',', skip_header = 1)
    if "cyt" in OptoChNames[i]:
        OptoChCyt = np.vstack((OptoChCyt, temp[:,2]))
    elif "jun" in OptoChNames[i]:
        OptoChJun = np.vstack((OptoChJun, temp[:,2]))
    else:
        print("variable mislabeled")

c450vC2Names = listdir(c450vC2Folder)

c450vC2Cyt = np.zeros((0,60))
c450vC2Jun = np.zeros((0,60))
for i in range(len(c450vC2Names)):
    source = join(c450vC2Folder, c450vC2Names[i])
    temp = np.genfromtxt(source, delimiter = ',', skip_header = 1)
    if "cyt" in c450vC2Names[i]:
        c450vC2Cyt = np.vstack((c450vC2Cyt, temp[:,2]))
    elif "jun" in c450vC2Names[i]:
        c450vC2Jun = np.vstack((c450vC2Jun, temp[:,2]))
    else:
        print("variable mislabeled")

c450vC5Names = listdir(c450vC5Folder)

c450vC5Cyt = np.zeros((0,60))
c450vC5Jun = np.zeros((0,60))
for i in range(len(c450vC5Names)):
    source = join(c450vC5Folder, c450vC5Names[i])
    temp = np.genfromtxt(source, delimiter = ',', skip_header = 1)
    if "cyt" in c450vC5Names[i]:
        c450vC5Cyt = np.vstack((c450vC5Cyt, temp[:,2]))
    elif "jun" in c450vC5Names[i]:
        c450vC5Jun = np.vstack((c450vC5Jun, temp[:,2]))
    else:
        print("variable mislabeled")

EyecupNames = listdir(EyecupFolder)

EyecupCyt = np.zeros((0,60))
EyecupJun = np.zeros((0,60))
for i in range(len(EyecupNames)):
    source = join(EyecupFolder, EyecupNames[i])
    temp = np.genfromtxt(source, delimiter = ',', skip_header = 1)
    if "cyt" in EyecupNames[i]:
        EyecupCyt = np.vstack((EyecupCyt, temp[:,2]))
    elif "jun" in EyecupNames[i]:
        EyecupJun = np.vstack((EyecupJun, temp[:,2]))
    else:
        print("variable mislabeled")

#Alright! Now it's time to do some maths!
#Let's pick all samples!
OptoShRated = np.zeros_like(OptoShCyt)
OptoChRated = np.zeros_like(OptoChCyt)
c450vC2Rated = np.zeros_like(c450vC2Cyt)
for i in range(OptoShCyt.shape[0]):
    OptoShRated[i,:] = OptoShJun[i,:]/OptoShCyt[i,:]
    OptoChRated[i,:] = OptoChJun[i,:]/OptoChCyt[i,:]
    c450vC2Rated[i,:] = c450vC2Jun[i,:]/c450vC2Cyt[i,:]
    
#We made more samples for eyecup
EyecupRated = np.zeros_like(EyecupCyt)
for i in range(EyecupCyt.shape[0]):
    EyecupRated[i,:] = EyecupJun[i,:]/EyecupCyt[i,:]

#Now averaging per row
OptoShRatedmean = []
OptoChRatedmean = []
c450vC2Ratedmean = []
EyecupRatedmean = []
OptoShRatedsdup = []
OptoChRatedsdup = []
c450vC2Ratedsdup = []
EyecupRatedsdup= []
OptoShRatedsdlow = []
OptoChRatedsdlow = []
c450vC2Ratedsdlow = []
EyecupRatedsdlow = []
#We need to make time list
timelist = []
for i in range(60):
    OptoShRatedmean.append(np.mean(OptoShRated[:,i]))
    OptoChRatedmean.append(np.mean(OptoChRated[:,i]))
    c450vC2Ratedmean.append(np.mean(c450vC2Rated[:,i]))
    EyecupRatedmean.append(np.mean(EyecupRated[:,i]))
    OptoShRatedsdup.append(np.mean(OptoShRated[:,i]) + np.std(OptoShRated[:,i]))
    OptoChRatedsdup.append(np.mean(OptoChRated[:,i]) + np.std(OptoChRated[:,i]))
    c450vC2Ratedsdup.append(np.mean(c450vC2Rated[:,i]) + np.std(c450vC2Rated[:,i]))
    EyecupRatedsdup.append(np.mean(EyecupRated[:,i]) + np.std(EyecupRated[:,i]))
    OptoShRatedsdlow.append(np.mean(OptoShRated[:,i]) - np.std(OptoShRated[:,i]))
    OptoChRatedsdlow.append(np.mean(OptoChRated[:,i]) - np.std(OptoChRated[:,i]))
    c450vC2Ratedsdlow.append(np.mean(c450vC2Rated[:,i]) - np.std(c450vC2Rated[:,i]))
    EyecupRatedsdlow.append(np.mean(EyecupRated[:,i]) - np.std(EyecupRated[:,i]))
    timelist.append(i*5.16)


import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['xtick.labelsize']=45
plt.rcParams['ytick.labelsize']=35
plt.rcParams.update({'font.size': 45})

#'Saving folder:
savefolder = r'C:\Users\ara\Documents\Thesis\OptoShroom3 Paper\RawData and Python Scripts\Translocation analysis'

fig,ax = plt.subplots(figsize = (9, 7), dpi=200)
ax.axvspan(60, 120, color = '#a4c3f6', alpha = 0.3)
#set_title("Translocation evolution (n=18) ratio junction/cytoplasm signal")
ax.plot(timelist, c450vC2Ratedmean, label = "OptoShroom3 c450v mutant 2", color = "#60D394",  linewidth=2.0)
ax.fill_between(timelist, c450vC2Ratedsdup,c450vC2Ratedsdlow , color='#60D394', alpha=.1)
ax.plot(timelist, OptoChRatedmean, label = "OptoShroom3" + "$\Delta$SD2", color = "#4f81c2",  linewidth=2.0)
ax.fill_between(timelist, OptoChRatedsdup,OptoChRatedsdlow , color='#4f81c2', alpha=.1)
ax.plot(timelist, EyecupRatedmean, label = "OptoShroom3 optic vesicles", color = "#820263",  linewidth=2.0)
ax.fill_between(timelist, EyecupRatedsdup,EyecupRatedsdlow , color='#820263', alpha=.1)
ax.plot(timelist, OptoShRatedmean, label = "OptoShroom3", color = "#dc756d",  linewidth=2.0)
ax.fill_between(timelist, OptoShRatedsdup,OptoShRatedsdlow , color='#dc756d', alpha=.1)
ax.set_ylabel("Signal ratio (junction/cytoplasm)")
ax.set_xlabel("time (seconds)")
ax.legend(fontsize =20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig(savefolder + r"//translocationRatioGraph.png", bbox_inches='tight')
plt.show()



#if we average them to their original val:
OptoShRatedmeanrel = []
OptoChRatedmeanrel = []
c450vC2Ratedmeanrel = []
EyecupRatedmeanrel=[]
for i in range(60):
    OptoShRatedmeanrel.append(np.mean(OptoShRatedmean[i]/OptoShRatedmean[0]))
    OptoChRatedmeanrel.append(np.mean(OptoChRatedmean[i]/OptoChRatedmean[0]))
    c450vC2Ratedmeanrel.append(np.mean(c450vC2Ratedmean[i]/c450vC2Ratedmean[0]))
    EyecupRatedmeanrel.append(np.mean(EyecupRatedmean[i]/EyecupRatedmean[0])) 

#Graph:
plt.figure(figsize = (9,7))
plt.title("Translocation evolution (n=18), ratio normalized ")
plt.plot(timelist, OptoShRatedmeanrel, label = "OptoShroom3", color = "blue")
plt.plot(timelist, OptoChRatedmeanrel, label = "OptomCherry", color = "red")
plt.plot(timelist, c450vC2Ratedmeanrel, label = "OptoShroom3 c450v mutant 2", color = "green")
plt.plot(timelist, EyecupRatedmeanrel, label = "OptoShroom3 Eyecups", color = "purple")
plt.xlabel("time (seconds)")
plt.ylabel("Signal ratio (junction/cytoplasm), relative to t = 0")
plt.axvspan(60, 120, color = "blue", alpha = 0.1, label = "stimulation period")
plt.legend(fontsize =20)
plt.savefig(savefolder + r"//translocationRatioNormalizedGraph.png", bbox_inches='tight')
plt.show()


#We need the same thing but with normalization to the initial value!

OptoShmeanJun = []
OptoChmeanJun = []
c450vC2meanJun = []
EyecupmeanJun = []
OptoShsdupJun = []
OptoChsdupJun = []
c450vC2sdupJun = []
EyecupsdupJun = []
OptoShsdlowJun = []
OptoChsdlowJun = []
c450vC2sdlowJun = []
EyecupsdlowJun = []
for i in range(60):
    OptoShmeanJun.append(np.mean(OptoShJun[:,i]/OptoShJun[:,0]))
    OptoChmeanJun.append(np.mean(OptoChJun[:,i]/OptoChJun[:,0]))
    c450vC2meanJun.append(np.mean(c450vC2Jun[:,i]/c450vC2Jun[:,0]))
    EyecupmeanJun.append(np.mean(EyecupJun[:,i]/EyecupJun[:,0]))
    OptoShsdupJun.append(np.mean(OptoShJun[:,i]/OptoShJun[:,0]) + np.std(OptoShJun[:,i]/OptoShJun[:,0]))
    OptoChsdupJun.append(np.mean(OptoChJun[:,i]/OptoChJun[:,0]) + np.std(OptoChJun[:,i]/OptoChJun[:,0]) )
    c450vC2sdupJun.append(np.mean(c450vC2Jun[:,i]/c450vC2Jun[:,0]) + np.std(c450vC2Jun[:,i]/c450vC2Jun[:,0]) )
    EyecupsdupJun.append(np.mean(EyecupJun[:,i]/EyecupJun[:,0]) + np.std(EyecupJun[:,i]/EyecupJun[:,0]))
    OptoShsdlowJun.append(np.mean(OptoShJun[:,i]/OptoShJun[:,0]) - np.std(OptoShJun[:,i]/OptoShJun[:,0]))
    OptoChsdlowJun.append(np.mean(OptoChJun[:,i]/OptoChJun[:,0]) - np.std(OptoChJun[:,i]/OptoChJun[:,0]) )
    c450vC2sdlowJun.append(np.mean(c450vC2Jun[:,i]/c450vC2Jun[:,0]) - np.std(c450vC2Jun[:,i]/c450vC2Jun[:,0]) )
    EyecupsdlowJun.append(np.mean(EyecupJun[:,i]/EyecupJun[:,0]) - np.std(EyecupJun[:,i]/EyecupJun[:,0]))


plt.figure(figsize = (9,7))
plt.plot(timelist, OptoShmeanJun, label = "OptoShroom3", color = "blue")
plt.plot(timelist, OptoChmeanJun, label = "OptomCherry", color = "red")
plt.plot(timelist, c450vC2meanJun, label = "OptoShroom3 c450v mutant 2", color = "green")
plt.plot(timelist, EyecupmeanJun, label = "OptoShroom3 Eyecups", color = "purple")
plt.axvspan(60, 120, color = "blue", alpha = 0.1)
plt.legend(fontsize =20)
plt.show()

fig,ax = plt.subplots(figsize = (9, 7), dpi=200)
ax.axvspan(60, 120, color = '#a4c3f6', alpha = 0.3, label = "stimulation period")
#set_title("Translocation evolution (n=18) ratio junction/cytoplasm signal")
ax.plot(timelist, c450vC2meanJun, label = "OptoShroom3 c450v mutant 2", color = "#60D394",  linewidth=2.0)
ax.fill_between(timelist, c450vC2sdupJun,c450vC2sdlowJun , color='#60D394', alpha=.1)
ax.plot(timelist, OptoChmeanJun, label = "OptoShroom3" + "$\Delta$SD2", color = "#4f81c2",  linewidth=2.0)
ax.fill_between(timelist, OptoChsdupJun,OptoChsdlowJun , color='#4f81c2', alpha=.1)
ax.plot(timelist, EyecupmeanJun, label = "OptoShroom3 optic vesicles", color = "#820263",  linewidth=2.0)
ax.fill_between(timelist, EyecupsdupJun,EyecupsdlowJun , color='#820263', alpha=.1)
ax.plot(timelist, OptoShmeanJun, label = "OptoShroom3", color = "#dc756d",  linewidth=2.0)
ax.fill_between(timelist, OptoShsdupJun,OptoShsdlowJun , color='#dc756d', alpha=.1)
ax.set_ylabel("Junction signal relative to t=0")
ax.set_xlabel("time (seconds)")
ax.legend(fontsize =20)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig(savefolder + r"//translocationRatioGraph.png", bbox_inches='tight')
plt.show()

#Same for cyt:

OptoShmeanCyt = []
OptoChmeanCyt = []
c450vC2meanCyt = []
EyecupmeanCyt = []
for i in range(60):
    OptoShmeanCyt.append(np.mean(OptoShCyt[:,i]/OptoShCyt[:,0]))
    OptoChmeanCyt.append(np.mean(OptoChCyt[:,i]/OptoChCyt[:,0]))
    c450vC2meanCyt.append(np.mean(c450vC2Cyt[:,i]/c450vC2Cyt[:,0]))
    EyecupmeanCyt.append(np.mean(EyecupCyt[:,i]/EyecupCyt[:,0]))

plt.figure(figsize = (9,7))
plt.title("Cytoplasm signal averages")
plt.plot(timelist, OptoShmeanCyt, label = "OptoShroom3", color = "blue")
plt.plot(timelist, OptoChmeanCyt, label = "OptomCherry", color = "red")
plt.plot(timelist, c450vC2meanCyt, label = "OptoShroom3 c450v mutant 2", color = "green")
plt.plot(timelist, EyecupmeanCyt, label = "OptoShroom3 Eyecups", color = "purple")
plt.xlabel("time (seconds)")
plt.ylabel("Cytoplasm signal, relative to t = 0")
plt.axvspan(60, 120, color = "blue", alpha = 0.1)
plt.legend(fontsize =20)
plt.show()


#Bleach correction:
OptoShmeanCorrected = []
OptoChmeanCorrected = []
c450vC2meanCorrected = []
EyecupmeanCorrected = []
for i in range(60):
    OptoShmeanCorrected.append(OptoShmeanJun[i]/OptoShmeanCyt[i])
    OptoChmeanCorrected.append(OptoChmeanJun[i]/OptoChmeanCyt[i])
    c450vC2meanCorrected.append(c450vC2meanJun[i]/c450vC2meanCyt[i])
    EyecupmeanCorrected.append(EyecupmeanJun[i]/EyecupmeanCyt[i])

plt.figure(figsize = (9,7))
plt.title("Averaged and corrected translocation evolution")
plt.plot(timelist, OptoShmeanCorrected, label = "OptoShroom3", color = "blue")
plt.plot(timelist, OptoChmeanCorrected, label = "OptomCherry", color = "red")
plt.plot(timelist, c450vC2meanCorrected, label = "OptoShroom3 c450v mutant", color = "green")
plt.plot(timelist, EyecupmeanCorrected, label = "OptoShroom3 Eyecups", color = "purple")
plt.axvspan(60, 120, color = "blue", alpha = 0.1)
plt.legend(fontsize =20)
plt.show()

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams.update({'font.size': 7})
plt.rcParams['pdf.fonttype'] = 42
#Finally, the graphs together that will make a supplementary figure for the paper.
fig,([ax1, ax2], [ax3, ax4]) = plt.subplots(2,2,figsize = (8, 7))
ax1.axvspan(60, 120, color = '#a4c3f6', alpha = 0.3, fill = "True", lw= 0)
#set_title("Translocation evolution (n=18) ratio junction/cytoplasm signal")
ax1.plot(timelist, OptoShmeanJun, label = "OptoShroom3 (MDCK)", color = "#dc756d",  linewidth=1)
ax1.fill_between(timelist, OptoShsdupJun, OptoShsdlowJun, color='#dc756d', alpha=.1)
ax1.plot(timelist, c450vC2meanJun, label = "OptoShroom3 C450V mutant (MDCK)", color = "#60D394",  linewidth=1)
ax1.fill_between(timelist, c450vC2sdupJun,c450vC2sdlowJun , color='#60D394', alpha=.1)
ax1.plot(timelist, OptoChmeanJun, label = "GFP-NShroom3-iLID + SspB-mCherry (MDCK)", color = "#4f81c2",  linewidth=1)
ax1.fill_between(timelist, OptoChsdupJun,OptoChsdlowJun , color='#4f81c2', alpha=.1)
ax1.set_ylabel("Junction signal relative to t=0")
ax1.set_xlabel("Time (seconds)")
ax1.legend(fontsize =6)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.axvspan(60, 120, color = '#a4c3f6', alpha = 0.3, fill = "True", lw= 0)
#set_title("Translocation evolution (n=18) ratio junction/cytoplasm signal")
ax2.plot(timelist, OptoShRatedmean, label = "OptoShroom3 (MDCK)", color = "#dc756d",  linewidth=1)
ax2.fill_between(timelist, OptoShRatedsdup,OptoShRatedsdlow , color='#dc756d', alpha=.1)
ax2.plot(timelist, c450vC2Ratedmean, label = "OptoShroom3 C450V mutant (MDCK)", color = "#60D394",  linewidth=1)
ax2.fill_between(timelist, c450vC2Ratedsdup,c450vC2Ratedsdlow , color='#60D394', alpha=.1)
ax2.plot(timelist, OptoChRatedmean, label = "GFP-NShroom3-iLID + SspB-mCherry (MDCK)", color = "#4f81c2",  linewidth=1)
ax2.fill_between(timelist, OptoChRatedsdup,OptoChRatedsdlow , color='#4f81c2', alpha=.1)
ax2.set_ylabel("Signal ratio (junction/cytoplasm)")
ax2.set_xlabel("Time (seconds)")
ax2.legend(fontsize =6)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3.axvspan(60, 120, color = '#a4c3f6', alpha = 0.3, fill = "True", lw= 0)
ax3.plot(timelist, EyecupmeanJun, label = "OptoShroom3 (Optic vesicle)", color = "#820263",  linewidth=1)
ax3.fill_between(timelist, EyecupsdupJun,EyecupsdlowJun , color='#820263', alpha=.1)
ax3.set_ylabel("Junction signal relative to t=0")
ax3.set_xlabel("Time (seconds)")
ax3.legend(fontsize =6)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax4.axvspan(60, 120, color = '#a4c3f6', alpha = 0.3, fill = "True", lw= 0)
ax4.plot(timelist, EyecupRatedmean, label = "OptoShroom3 (Optic vesicle)", color = "#820263",  linewidth=1)
ax4.fill_between(timelist, EyecupRatedsdup,EyecupRatedsdlow , color='#820263', alpha=.1)
ax4.set_ylabel("Signal ratio (junction/cytoplasm)")
ax4.set_xlabel("Time (seconds)")
ax4.legend(fontsize =6)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
fig.savefig(savefolder + r"//translocationSuppfigure.pdf", bbox_inches='tight', transparent=True)
plt.show()

len(c450vC2Jun)
len(EyecupJun)


#Let's export all of this data by panel:
#a
import pandas as pd
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"
dataexp = pd.DataFrame(list(zip(timelist, OptoShmeanJun, OptoShsdupJun, OptoShsdlowJun, c450vC2meanJun, c450vC2sdupJun, OptoShsdlowJun, OptoChmeanJun, OptoChsdupJun,OptoChsdlowJun)), 
                       columns = ["Time (min)", "OptoShroom3 average", "OptoShroom3 average + standard deviation", "OptoShroom3 average - standard deviation", "OptoShroom3 no CShroom3 average", "OptoShroom3 no CShroom3 average + standard deviation", "OptoShroom3 no CShroom3 average - standard deviation", "OptoShroom3 C450V average", "OptoShroom3 C450V average + standard deviation", "OptoShroom3 C450V average - standard deviation"])

dataexp.to_csv(filepath + '\\supfig2a.csv', index = False)
#b
dataexp = pd.DataFrame(list(zip(timelist, OptoShRatedmean, OptoShRatedsdup, OptoShRatedsdlow, c450vC2Ratedmean, c450vC2Ratedsdup, c450vC2Ratedsdlow, OptoChRatedmean, OptoChRatedsdup, OptoChRatedsdlow)), 
                       columns = ["Time (min)", "OptoShroom3 average", "OptoShroom3 average + standard deviation", "OptoShroom3 average - standard deviation", "OptoShroom3 no CShroom3 average", "OptoShroom3 no CShroom3 average + standard deviation", "OptoShroom3 no CShroom3 average - standard deviation", "OptoShroom3 C450V average", "OptoShroom3 C450V average + standard deviation", "OptoShroom3 C450V average - standard deviation"])

dataexp.to_csv(filepath + '\\supfig2b.csv', index = False)
#c
dataexp = pd.DataFrame(list(zip(timelist, EyecupmeanJun, EyecupsdupJun, EyecupsdlowJun)), 
                       columns = ["Time (min)", "OptoShroom3 optic vesicle average", "OptoShroom3 optic vesicle average + standard deviation", "OptoShroom3 optic vesicle average - standard deviation"])

dataexp.to_csv(filepath + '\\supfig2c.csv', index = False)
#d
dataexp = pd.DataFrame(list(zip(timelist, EyecupRatedmean, EyecupRatedsdup, EyecupRatedsdlow)), 
                       columns = ["Time (min)", "OptoShroom3 optic vesicle average", "OptoShroom3 optic vesicle average + standard deviation", "OptoShroom3 optic vesicle average - standard deviation"])

dataexp.to_csv(filepath + '\\supfig2d.csv', index = False)