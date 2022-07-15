# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:15:32 2022

@author: ara
"""

# Some quick code to check if classification by cell layer helps us understand constriction dynamics

#This first trial will be with apical signal only

#I'm going to get the originally segmented file, polish it and get the top apical projection
#Then on the first frame before stimulation I will classify the cells:
    #Stimulated: more than 50% of apical area in stimulation area
    #Layer1: Cells in contact with stimulated
    #Layer2: Cells in contact with Layer 1
    #Layer3: Cells in contact with Layer 2
    #Layer4: Cells in contact with Layer 3
#I dont think there will be any other useful layers;

#Trial for one sample only, to be automated with the others:
import numpy as np
import matplotlib.pyplot as plt
from tifffile import imwrite, TiffFile
from os.path import join
import scipy.ndimage as ndi
from skimage import filters, measure, morphology, transform
from scipy import stats
import math
import matplotlib
    

path = r"C:\Users\ara\Documents\Python Scripts\3D segmentation MDCK cells\Reversible Stimulation Tests\CellLayerClasification tests\Samples"
graphpath = r'C:\Users\ara\Documents\Python Scripts\3D segmentation MDCK cells\Reversible Stimulation Tests\CellLayerClasification tests\Graphs'
xyres = 0.2071602
zres = 0.7

#Final data I will put here:

    
#colormap we will use for display: 
cmap = np.random.rand ( 256,3)
cmap[0,:] = 1
color = matplotlib.colors.ListedColormap(cmap)

from os import listdir
filenames = listdir(path)

finalstimulated = []
finallayer1 = []
finallayer2 = []
finallayer3 = []
for n in filenames:
    print("We are on sample " + n)
    samplename = n
    
    samplePath = join(path, samplename)
    
    with TiffFile(samplePath) as tif:
        sample = tif.asarray()
    
    frames, z, chan, y, x = sample.shape
    
    cells3d = sample[:,:,1,:,:]
    
    
    cells3d = cells3d-np.amin(cells3d)
    
    cells3dSmooth = np.zeros_like(cells3d)
    borders = np.zeros_like(cells3d[0,:,:,:]) #(t,z,channel1,x,y) 
    borders[0,:,:] = 1
    borders[:,0,:] = 1
    borders[:,:,0] = 1
    borders[-1,:,:] = 1
    borders[:,-1,:] = 1
    borders[:,:,-1] = 1
    
    #The smooth file is all we need to get the apical cuts
    #Then, we will get the selected list of the cells we are looking into by checking only the cells that are there all the time. 
    print("Smoothing")
    skipcells = []
    for t in range(frames):
        print("Timepoint = " + str(t))
        for i in np.unique(cells3d[t,:,:,:]):
            if i == 0:
                continue
            #print("We are on cell: " + str(i))
            currentcell = cells3d[t,:,:,:] == i
            currentcell = ndi.binary_erosion(currentcell, iterations = 2)
            currentcell = ndi.binary_dilation(currentcell, iterations = 2)
            cells3dSmooth[t,:,:,:] = np.where(currentcell == 1 , i, cells3dSmooth[t,:,:,:])
        hits = borders*cells3dSmooth[t,:,:,:]
        for b in np.unique(hits):
            if b in skipcells:
                pass
            else:
                skipcells.append(b)
        #Now check all the cells
    
        
    allcells = np.unique(cells3d)
    
    for i in allcells:
        print(i)
        for t in range(frames):
            if i not in cells3dSmooth[t,:,:,:]:
                skipcells.append(i)
                break
            else:
                continue
    
    #Now the list of cells to look into:
    cellskept = []
    for i in allcells:
        if i not in skipcells:
            cellskept.append(i)
    
    #Ok now just generate the apical image!
    print("projecting")
    apicalcuts = np.zeros((frames,x,y))
    
    for t in range(frames):
        print(t)
        for a in range(y):
            for b in range(x):
                for c in range(z):
                    if cells3dSmooth[t,z-1-c,a,b] == 0:
                        pass
                    else:
                        apicalcuts[t,a,b] = cells3dSmooth[t,z-1-c,a,b]
                        break
    
    #Ok and now on to the clasification, we are playing on timepoint 4
    
    
    stimulated = []
    stimulatedall = []
    #stimulated cells will be only those whose 50% area or more is within the area
    stimcandidates = list(np.unique(apicalcuts[3,192:320,192:320]))
    stimcandidates.pop(0)
    for i in stimcandidates: 
        currentcell = apicalcuts[3,:,:] == i
        partinside = apicalcuts[3,192:320,192:320] == i
        currentcellarea = np.sum(currentcell)
        partinsidearea = np.sum(partinside)
        ratio = partinsidearea/currentcellarea
        if ratio>0.5:
            stimulatedall.append(i)
            if i in cellskept:
                stimulated.append(i)
            
    #Now on to the rest of the layers!
    
    layer1 = []
    layer2 = []
    layer3 = []
    layer4 = []
    #layernall versions are to make the images!
    layer1all = []
    layer2all = []
    layer3all = []
    layer4all = []
    
    #1 Get the cells within the first layer
    for i in stimulatedall:
        currentcell = apicalcuts[3,:,:] == i
        currentcelldil = ndi.binary_dilation(currentcell)
        perimeter = currentcelldil^currentcell
        neighbours = perimeter*apicalcuts
        candidates = np.unique(neighbours)
        for a in candidates:
            if a == 0: 
                continue
            elif a in stimulatedall:
                continue
            elif a in layer1all:
                continue
            else:
                times = np.count_nonzero(neighbours == a)
                if times > 9:
                    #reliable nighbour
                    layer1all.append(a)
                    if a in cellskept:
                        layer1.append(a)
        
    #Cool! Now on to layer 2:
    
    for i in layer1all:
        currentcell = apicalcuts[3,:,:] == i
        currentcelldil = ndi.binary_dilation(currentcell)
        perimeter = currentcelldil^currentcell
        neighbours = perimeter*apicalcuts
        candidates = np.unique(neighbours)
        for a in candidates:
            if a == 0: 
                continue
            elif a in stimulatedall:
                continue
            elif a in layer1all:
                continue
            elif a in layer2all:
                continue
            else:
                times = np.count_nonzero(neighbours == a)
                if times > 9:
                    #reliable nighbour
                    layer2all.append(a)
                    if a in cellskept:
                        layer2.append(a)
    
    
    #Cool! Now on to layer 3:
    
    for i in layer2all:
        currentcell = apicalcuts[3,:,:] == i
        currentcelldil = ndi.binary_dilation(currentcell)
        perimeter = currentcelldil^currentcell
        neighbours = perimeter*apicalcuts
        candidates = np.unique(neighbours)
        for a in candidates:
            if a == 0: 
                continue
            elif a in stimulatedall:
                continue
            elif a in layer1all:
                continue
            elif a in layer2all:
                continue
            elif a in layer3all:
                continue
            else:
                times = np.count_nonzero(neighbours == a)
                if times > 9:
                    #reliable nighbour
                    layer3all.append(a)
                    if a in cellskept:
                        layer3.append(a)
    
            
    #No need for layer 4 really
    
    
    #Cooool! Let's make some images with the actual layers?
    
    
    stimulatedim = np.zeros_like(apicalcuts[3,:,:])
    for i in stimulatedall:
        currentcell = apicalcuts[3,:,:]==i
        stimulatedim = np.where(currentcell==True, i, stimulatedim)
    
    layer1im = np.zeros_like(apicalcuts[3,:,:])
    for i in layer1all:
        currentcell = apicalcuts[3,:,:]==i
        layer1im = np.where(currentcell==True, i, layer1im)
    
    layer2im = np.zeros_like(apicalcuts[3,:,:])
    for i in layer2all:
        currentcell = apicalcuts[3,:,:]==i
        layer2im = np.where(currentcell==True, i, layer2im)
    
    layer3im = np.zeros_like(apicalcuts[3,:,:])
    for i in layer3all:
        currentcell = apicalcuts[3,:,:]==i
        layer3im = np.where(currentcell==True, i, layer3im)
        
    

    plt.figure(figsize=(8,8))
    plt.suptitle(samplename)
    plt.subplot(2,2,1)
    plt.title("Stimulated cells")
    plt.imshow(stimulatedim, cmap=color, interpolation = "none")
    plt.subplot(2,2,2)
    plt.title("Layer1")
    plt.imshow(layer1im, cmap=color, interpolation = "none")
    plt.subplot(2,2,3)
    plt.title("layer2")
    plt.imshow(layer2im, cmap=color, interpolation = "none")
    plt.subplot(2,2,4)
    plt.title("layer3")
    plt.imshow(layer3im, cmap=color, interpolation = "none")
    plt.savefig(graphpath +r"\\" + samplename[:-4] + r"layerimages.png", bbox_inches='tight')
    plt.show()
    
    #Beautiful!!! Can we see some areas?
    stimareas = []
    for i in stimulated:
        temp = []
        for t in range(frames):
            currentcell = apicalcuts[t,:,:] == i
            temp.append(np.sum(currentcell)*xyres*xyres)
        stimareas.append(temp)
    
    layer1areas = []
    for i in layer1:
        temp = []
        for t in range(frames):
            currentcell = apicalcuts[t,:,:] == i
            temp.append(np.sum(currentcell)*xyres*xyres)
        layer1areas.append(temp)
        
    layer2areas = []
    for i in layer2:
        temp = []
        for t in range(frames):
            currentcell = apicalcuts[t,:,:] == i
            temp.append(np.sum(currentcell)*xyres*xyres)
        layer2areas.append(temp)
    
    layer3areas = []
    for i in layer3:
        temp = []
        for t in range(frames):
            currentcell = apicalcuts[t,:,:] == i
            temp.append(np.sum(currentcell)*xyres*xyres)
        layer3areas.append(temp)
        
        
    #And now averages!!!
    stimave = []
    layer1ave = []
    layer2ave = []
    layer3ave = []
    for t in range(frames):
        print(t)
        stimtemp = []
        for i in range(len(stimareas)):
            stimtemp.append(stimareas[i][t])
        stimave.append(np.mean(stimtemp))
        layer1temp = []
        for i in range(len(layer1areas)):
            layer1temp.append(layer1areas[i][t])
        layer1ave.append(np.mean(layer1temp))
        layer2temp = []
        for i in range(len(layer2areas)):
            layer2temp.append(layer2areas[i][t])
        layer2ave.append(np.mean(layer2temp))
        layer3temp = []
        for i in range(len(layer3areas)):
            layer3temp.append(layer3areas[i][t])
        layer3ave.append(np.mean(layer3temp))
        
    timepoints = [i*10 for i in range(frames)]
    
    plt.figure(figsize = (7,5))
    plt.plot(timepoints, stimave, label = "stimulated")
    plt.plot(timepoints, layer1ave, label = "layer1")
    plt.plot(timepoints, layer2ave, label = "layer2")
    plt.plot(timepoints, layer3ave, label = "layer3")
    plt.legend()
    plt.xlabel("time (min)")
    plt.ylabel("Apical area "+ r"$(μm^{2})$")
    plt.savefig(graphpath +r"\\" + samplename[:-4] + r"values.png", bbox_inches='tight')
    plt.show()
    
    #Finally we just need to output all the apical measurements:
    
    finalstimulated = finalstimulated + stimareas
    finallayer1 = finallayer1 + layer1areas
    finallayer2 = finallayer2 + layer2areas
    finallayer3 = finallayer3 + layer3areas

#That should be good to go!
#Layer by layer analysis is ready to be executed!
    
#Time to get averages 
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

stimlist = getAverage(finalstimulated)
layer1list = getAverage(finallayer1)
layer2list = getAverage(finallayer2)
layer3list = getAverage(finallayer3)

#Getting N for figure foot!
len(finalstimulated)
len(finallayer1)
len(finallayer2)
len(finallayer3)


#Now the plot!

plt.figure(figsize = (7,5))
plt.plot(timepoints, stimlist[0], label = "stimulated")
plt.plot(timepoints, layer1list[0], label = "layer1")
plt.plot(timepoints, layer2list[0], label = "layer2")
plt.plot(timepoints, layer3list[0], label = "layer3")
plt.legend()
plt.xlabel("time (min)")
plt.ylabel("Apical area "+ r"$(μm^{2})$")
plt.savefig(graphpath +r"\\"  + r"layerclasification.png", bbox_inches='tight')
plt.show()

#Now the relative graph!



relstimfinal = [[finalstimulated[i][n]/finalstimulated[i][3] for n in range(len(finalstimulated[i]))] for i in range(len(finalstimulated))]
rellayer1final = [[finallayer1[i][n]/finallayer1[i][3] for n in range(len(finallayer1[i]))] for i in range(len(finallayer1))]
rellayer2final = [[finallayer2[i][n]/finallayer2[i][3] for n in range(len(finallayer2[i]))] for i in range(len(finallayer2))]
rellayer3final = [[finallayer3[i][n]/finallayer3[i][3] for n in range(len(finallayer3[i]))] for i in range(len(finallayer3))]

relstimlist = getAverage(relstimfinal)
rellayer1list = getAverage(rellayer1final)
rellayer2list = getAverage(rellayer2final)
rellayer3list = getAverage(rellayer3final)

stimstart = 30
stimend2h = 150

plt.figure(figsize = (7,5))
plt.axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
plt.plot(timepoints, relstimlist[0], label = "stimulated", color = '#dc756d')
plt.fill_between(timepoints, relstimlist[1],relstimlist[2], color = '#dc756d', alpha= 0.05)
plt.plot(timepoints, rellayer1list[0], label = "layer 1" , color = '#4f81c2')
plt.fill_between(timepoints, rellayer1list[1], rellayer1list[2], color = '#4f81c2', alpha= 0.05)
plt.plot(timepoints, rellayer2list[0], label = "layer 2" , color = '#820263')
plt.fill_between(timepoints, rellayer2list[1], rellayer2list[2], color = '#820263', alpha= 0.05)
plt.plot(timepoints, rellayer3list[0], label = "layer 3" , color = '#60D394')
plt.fill_between(timepoints, rellayer3list[1], rellayer3list[2], color = '#60D394', alpha= 0.05)
#plt.ylim(0,1.2)
plt.legend()
plt.xlabel("time (min)")
plt.ylabel("Relative apical area ")
plt.savefig(graphpath +r"\\"  + r"rellayerclasification.png", bbox_inches='tight')
plt.show()

#Ok and the figure for the paper: 

#We have to make the image for the 4 layers, I decided to do it with sample 202200203_2hstim3dseg

samplename = filenames[3]

#We just remake the loop to get all data again:
samplePath = join(path, samplename)

with TiffFile(samplePath) as tif:
    sample = tif.asarray()

frames, z, chan, y, x = sample.shape

cells3d = sample[:,:,1,:,:]


cells3d = cells3d-np.amin(cells3d)

cells3dSmooth = np.zeros_like(cells3d)
borders = np.zeros_like(cells3d[0,:,:,:]) #(t,z,channel1,x,y) 
borders[0,:,:] = 1
borders[:,0,:] = 1
borders[:,:,0] = 1
borders[-1,:,:] = 1
borders[:,-1,:] = 1
borders[:,:,-1] = 1

#The smooth file is all we need to get the apical cuts
#Then, we will get the selected list of the cells we are looking into by checking only the cells that are there all the time. 
print("Smoothing")
skipcells = []
for t in range(frames):
    print("Timepoint = " + str(t))
    for i in np.unique(cells3d[t,:,:,:]):
        if i == 0:
            continue
        #print("We are on cell: " + str(i))
        currentcell = cells3d[t,:,:,:] == i
        currentcell = ndi.binary_erosion(currentcell, iterations = 2)
        currentcell = ndi.binary_dilation(currentcell, iterations = 2)
        cells3dSmooth[t,:,:,:] = np.where(currentcell == 1 , i, cells3dSmooth[t,:,:,:])
    hits = borders*cells3dSmooth[t,:,:,:]
    for b in np.unique(hits):
        if b in skipcells:
            pass
        else:
            skipcells.append(b)
    #Now check all the cells

    
allcells = np.unique(cells3d)

for i in allcells:
    print(i)
    for t in range(frames):
        if i not in cells3dSmooth[t,:,:,:]:
            skipcells.append(i)
            break
        else:
            continue

#Now the list of cells to look into:
cellskept = []
for i in allcells:
    if i not in skipcells:
        cellskept.append(i)

#Ok now just generate the apical image!
print("projecting")
apicalcuts = np.zeros((frames,x,y))

for t in range(frames):
    print(t)
    for a in range(y):
        for b in range(x):
            for c in range(z):
                if cells3dSmooth[t,z-1-c,a,b] == 0:
                    pass
                else:
                    apicalcuts[t,a,b] = cells3dSmooth[t,z-1-c,a,b]
                    break

#Ok and now on to the clasification, we are playing on timepoint 4


stimulated = []
stimulatedall = []
#stimulated cells will be only those whose 50% area or more is within the area
stimcandidates = list(np.unique(apicalcuts[3,192:320,192:320]))
stimcandidates.pop(0)
for i in stimcandidates: 
    currentcell = apicalcuts[3,:,:] == i
    partinside = apicalcuts[3,192:320,192:320] == i
    currentcellarea = np.sum(currentcell)
    partinsidearea = np.sum(partinside)
    ratio = partinsidearea/currentcellarea
    if ratio>0.5:
        stimulatedall.append(i)
        if i in cellskept:
            stimulated.append(i)
        
#Now on to the rest of the layers!

layer1 = []
layer2 = []
layer3 = []
layer4 = []
#layernall versions are to make the images!
layer1all = []
layer2all = []
layer3all = []
layer4all = []

#1 Get the cells within the first layer
for i in stimulatedall:
    currentcell = apicalcuts[3,:,:] == i
    currentcelldil = ndi.binary_dilation(currentcell)
    perimeter = currentcelldil^currentcell
    neighbours = perimeter*apicalcuts
    candidates = np.unique(neighbours)
    for a in candidates:
        if a == 0: 
            continue
        elif a in stimulatedall:
            continue
        elif a in layer1all:
            continue
        else:
            times = np.count_nonzero(neighbours == a)
            if times > 9:
                #reliable nighbour
                layer1all.append(a)
                if a in cellskept:
                    layer1.append(a)
    
#Cool! Now on to layer 2:

for i in layer1all:
    currentcell = apicalcuts[3,:,:] == i
    currentcelldil = ndi.binary_dilation(currentcell)
    perimeter = currentcelldil^currentcell
    neighbours = perimeter*apicalcuts
    candidates = np.unique(neighbours)
    for a in candidates:
        if a == 0: 
            continue
        elif a in stimulatedall:
            continue
        elif a in layer1all:
            continue
        elif a in layer2all:
            continue
        else:
            times = np.count_nonzero(neighbours == a)
            if times > 9:
                #reliable nighbour
                layer2all.append(a)
                if a in cellskept:
                    layer2.append(a)


#Cool! Now on to layer 3:

for i in layer2all:
    currentcell = apicalcuts[3,:,:] == i
    currentcelldil = ndi.binary_dilation(currentcell)
    perimeter = currentcelldil^currentcell
    neighbours = perimeter*apicalcuts
    candidates = np.unique(neighbours)
    for a in candidates:
        if a == 0: 
            continue
        elif a in stimulatedall:
            continue
        elif a in layer1all:
            continue
        elif a in layer2all:
            continue
        elif a in layer3all:
            continue
        else:
            times = np.count_nonzero(neighbours == a)
            if times > 9:
                #reliable nighbour
                layer3all.append(a)
                if a in cellskept:
                    layer3.append(a)
    

#Here I will try to make the image either with just different colours.
#Or with monochrome colors each, but then I have to make a cutom color map?
stimulatedim = np.zeros_like(apicalcuts[3,:,:])
for i in stimulatedall:
    currentcell = apicalcuts[3,:,:]==i
    stimulatedim = np.where(currentcell==True, i, stimulatedim)

layer1im = np.zeros_like(apicalcuts[3,:,:])
for i in layer1all:
    currentcell = apicalcuts[3,:,:]==i
    layer1im = np.where(currentcell==True, i, layer1im)

layer2im = np.zeros_like(apicalcuts[3,:,:])
for i in layer2all:
    currentcell = apicalcuts[3,:,:]==i
    layer2im = np.where(currentcell==True, i, layer2im)

layer3im = np.zeros_like(apicalcuts[3,:,:])
for i in layer3all:
    currentcell = apicalcuts[3,:,:]==i
    layer3im = np.where(currentcell==True, i, layer3im)
    
#We need to make masks for layered image!
stimmask = np.ma.masked_where(stimulatedim==0, stimulatedim)
layer1mask = np.ma.masked_where(layer1im==0, layer1im)
layer2mask = np.ma.masked_where(layer2im==0, layer2im)
layer3mask = np.ma.masked_where(layer3im==0, layer3im)

#colormap for each layer!
#We need to work it out from our main colours, which are in hex, 
from PIL import ImageColor

#if we make it tidy, we can generate the 4 colormaps:
def makeColorMap(color, factor = 0.7):
    #From a HEX color:
    rgbcol = ImageColor.getcolor(color, "RGB")
    #scale it
    rgbscaled = [i/256 for i in rgbcol]
    rgbcolend = [i*factor for i in rgbscaled]
    cmap1 = np.array((np.linspace(rgbcolend[0],rgbscaled[0], 256), np.linspace(rgbcolend[1],rgbscaled[1], 256), np.linspace(rgbcolend[2],rgbscaled[2], 256)))
    #But this map has dimensions in the wrong order:
    cmap2 = np.moveaxis(cmap1, 0, -1)
    #readyfor matplotlib:
    cmap = matplotlib.colors.ListedColormap(cmap2)
    return cmap
       
    
stimcol = makeColorMap("#dc756d")
lay1col = makeColorMap('#4f81c2')
lay2col = makeColorMap('#820263')
lay3col = makeColorMap('#60D394')




plt.figure(figsize=(8,8), frameon=False)
plt.suptitle(samplename)
plt.imshow(stimmask, cmap=stimcol, interpolation = "none")
plt.imshow(layer1mask, cmap=lay1col, interpolation = "none")
plt.imshow(layer2mask, cmap=lay2col, interpolation = "none")
plt.imshow(layer3mask, cmap=lay3col, interpolation = "none")
plt.savefig(graphpath +r"\\" + samplename[:-4] + r"layered.png", bbox_inches='tight')
plt.show()



#Now for the final image!

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['xtick.labelsize']=45
plt.rcParams['ytick.labelsize']=35
plt.rcParams.update({'font.size': 45})
from matplotlib.lines import Line2D

fig, [ax1, ax2] = plt.subplots(1,2, figsize = (28,12), dpi=200)
ax1.imshow(stimmask, cmap=stimcol, interpolation = "none", label = "Stimulated")
ax1.imshow(layer1mask, cmap=lay1col, interpolation = "none", label = "Layer 1")
ax1.imshow(layer2mask, cmap=lay2col, interpolation = "none", label = "Layer 2")
ax1.imshow(layer3mask, cmap=lay3col, interpolation = "none", label = "Layer 3")
legend_elements = [Line2D([0], [0], color="#dc756d", lw=6, label='Stimulated'),
                   Line2D([0], [0], color='#4f81c2', lw=6, label='Layer 1'),
                   Line2D([0], [0], color='#820263', lw=6, label='Layer 2'),
                   Line2D([0], [0], color='#60D394', lw=6, label='Layer 3')]
ax1.legend(handles=legend_elements, loc = "upper right", fontsize = 32)
ax1.axis('off')
ax2.axvspan(stimstart, stimend2h, color='#a4c3f6', alpha = 0.2)
ax2.plot(timepoints, relstimlist[0], label = "Stimulated", color = '#dc756d', linewidth = 6)
ax2.fill_between(timepoints, relstimlist[1],relstimlist[2], color = '#dc756d', alpha= 0.1)
ax2.plot(timepoints, rellayer1list[0], label = "Layer 1" , color = '#4f81c2', linewidth = 6)
ax2.fill_between(timepoints, rellayer1list[1], rellayer1list[2], color = '#4f81c2', alpha= 0.1)
ax2.plot(timepoints, rellayer2list[0], label = "Layer 2" , color = '#820263', linewidth = 6)
ax2.fill_between(timepoints, rellayer2list[1], rellayer2list[2], color = '#820263', alpha= 0.1)
ax2.plot(timepoints, rellayer3list[0], label = "Layer 3" , color = '#60D394', linewidth = 6)
ax2.fill_between(timepoints, rellayer3list[1], rellayer3list[2], color = '#60D394', alpha= 0.1)
#plt.ylim(0,1.2)
ax2.legend(loc = "upper right", fontsize = 32)
ax2.set_xlabel("Time (min)")
ax2.set_ylabel("Relative apical area")
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
plt.savefig(graphpath +r"\\"  + r"layerclasfig2.png", bbox_inches='tight')
plt.show()

#Exporting the data:
import pandas as pd
filepath = r"\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\SourceData\Individual Data"
dataexp = pd.DataFrame(list(zip(timepoints, relstimlist[0], relstimlist[1], relstimlist[2], rellayer1list[0], rellayer1list[1], rellayer1list[2], rellayer2list[0], rellayer2list[1], rellayer2list[2], rellayer3list[0], rellayer3list[1], rellayer3list[2])), 
                       columns = ["Time (min)", "Area stimulated", "Area stimulated + standard deviation", "Area stimulated - standard deviation", "Area layer 1", "Area layer 1 + standard deviation", "Area layer 1 - standard deviation", "Area layer 2", "Area layer 2 + standard deviation", "Area layer 2 - standard deviation", "Area layer 3", "Area layer 3 + standard deviation", "Area layer 3 - standard deviation"])

dataexp.to_csv(filepath + '\\layeranalisisconstriction.csv', index = False)
