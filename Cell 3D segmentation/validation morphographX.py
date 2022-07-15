# -*- coding: utf-8 -*-
"""
Created on Mon May  2 15:28:00 2022

@author: ara
"""

###Code to analyze 3D segmented data from morphographX and my own 3D segmentation

#The idea is to get all the segmented cells in one timepoint. Exclude those touching the border,
#Find the same cell in the two images and measure volume, height, apical and basal area.

#It must be noted that the samples from Morphograph are not isometric. zresolution is 0.7 and xyres = 0.207...

#Also, they came flipped on y, so I unflipped them on fiji so I can trace cells from one to the other. Sample prep for this analysis was don manually.

#My samples, on the other hand are raw from segmentation, so they needed to have a smoothing I do to all and to eliminate those touching borders.


#Samples from morphograph will be labeled "morpho"

#Samples from my 3D pipeline will be labeled as "pipe"

#First step load the images
import numpy as np
import matplotlib.pyplot as plt
from tifffile import imwrite, TiffFile
from os.path import join
import scipy.ndimage as ndi
from skimage import filters, measure, morphology, transform
from scipy import stats
import math

path = r"H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\Validations"
xyres = 0.2071602
zres = 0.7

morphoEasyName = r"flippedEasy_t1.tif"
morphoHardName = r"flippedHard_t1.tif"
pipeEasyName = r"3dSegEasy.tif"
pipeHardName = r"3dSegHard.tif"

morphoEasyPath = join(path, morphoEasyName)

with TiffFile(morphoEasyPath) as tif:
    morphoEasy = tif.asarray()

morphoHardPath = join(path, morphoHardName)

with TiffFile(morphoHardPath) as tif:
    morphoHard = tif.asarray()


pipeEasyPath = join(path, pipeEasyName)

with TiffFile(pipeEasyPath) as tif:
    pipeEasy = tif.asarray()

pipeHardPath = join(path, pipeHardName)

with TiffFile(pipeHardPath) as tif:
    pipeHard = tif.asarray()
pipeHard.shape

#Cool! Now we clean up the pipeline data:
pipeEasy = pipeEasy-np.amin(pipeEasy)
pipeEasySmooth = np.zeros_like(pipeEasy)
pipeEasyClean = np.zeros_like(pipeEasy)
cells3d = np.copy(pipeEasy)
borders = np.zeros_like(cells3d) #(t,z,channel1,x,y) 
borders[0,:,:] = 1
borders[:,0,:] = 1
borders[:,:,0] = 1
borders[-1,:,:] = 1
borders[:,-1,:] = 1
borders[:,:,-1] = 1

for i in np.unique(pipeEasy):
    print("We are on cell: " + str(i))
    currentcell = cells3d == i
    currentcell = ndi.binary_erosion(currentcell, iterations = 2)
    currentcell = ndi.binary_dilation(currentcell, iterations = 2)
    pipeEasySmooth = np.where(currentcell == 1 , i, pipeEasySmooth)
hits = borders*pipeEasySmooth
skip = np.unique(hits)
for i in np.unique(pipeEasySmooth):   
    currentcell = pipeEasySmooth == i
    if i in skip:
        continue
    else:
        print("We are on cell: " + str(i))
        pipeEasyClean = np.where(currentcell == 1 , i, pipeEasyClean)


#Same for hard:
pipeHard = pipeHard-np.amin(pipeHard)
pipeHardSmooth = np.zeros_like(pipeHard)
pipeHardClean = np.zeros_like(pipeHard)
cells3d = np.copy(pipeHard)
borders = np.zeros_like(cells3d) #(t,z,channel1,x,y) 
borders[0,:,:] = 1
borders[:,0,:] = 1
borders[:,:,0] = 1
borders[-1,:,:] = 1
borders[:,-1,:] = 1
borders[:,:,-1] = 1

for i in np.unique(pipeHard):
    print("We are on cell: " + str(i))
    currentcell = cells3d == i
    currentcell = ndi.binary_erosion(currentcell, iterations = 2)
    currentcell = ndi.binary_dilation(currentcell, iterations = 2)
    pipeHardSmooth = np.where(currentcell == 1 , i, pipeHardSmooth)
hits = borders*pipeHardSmooth
skip = np.unique(hits)
for i in np.unique(pipeHardSmooth):   
    currentcell = pipeHardSmooth == i
    if i in skip:
        continue
    else:
        print("We are on cell: " + str(i))
        pipeHardClean = np.where(currentcell == 1 , i, pipeHardClean)

plt.imshow(morphoHard[12,:,:])
##Coool! Now we have to go cell by cell, track to the one in the same position (max overlap), and measure volume height and areas

#Before entering the loop we do the apical and basal area measurements:

z, x, y = pipeEasySmooth.shape
pipeEasyApical = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if pipeEasySmooth[z-1-c,a,b] == 0:
                pass
            else:
                pipeEasyApical[a,b] = pipeEasySmooth[z-1-c,a,b]
                break    
plt.imshow(pipeEasyApical)

#And basal:
pipeEasyBasal = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if pipeEasySmooth[c,a,b] == 0:
                pass
            else:
                pipeEasyBasal[a,b] = pipeEasySmooth[c,a,b]
                break    
plt.imshow(pipeEasyBasal) 

#Maybe we have to polish a bit morpho segmentation by eliminating thin slices that were unassigned
morphoEasySmooth = np.zeros_like(morphoEasy)
for i in np.unique(morphoEasy):
    if i == 0:
        continue
    currentcell = morphoEasy == i
    currentcell = currentcell.astype('uint8')
    Props = measure.regionprops(currentcell)
    vol = Props[0].area*xyres*xyres*zres
    if vol>400:
        morphoEasySmooth = np.where(currentcell==1, i, morphoEasySmooth)

#And now apicals and basals:
z, x, y = morphoEasySmooth.shape
morphoEasyApical = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if morphoEasySmooth[z-1-c,a,b] == 0:
                pass
            else:
                morphoEasyApical[a,b] = morphoEasySmooth[z-1-c,a,b]
                break    
plt.imshow(morphoEasyApical)

#And basal:
morphoEasyBasal = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if morphoEasySmooth[c,a,b] == 0:
                pass
            else:
                morphoEasyBasal[a,b] = morphoEasySmooth[c,a,b]
                break    
plt.imshow(morphoEasyBasal) 

pipevols = []
pipeapicals = []
pipebasals = []
pipelengths = []
morphovols = []
morphoapicals = []
morphobasals = []
morpholengths = []
for i in np.unique(pipeEasyClean):
    if i ==0:
        continue
    else:
        #First, find parallel cell:
        print("on cell " + str(i))
        pipecell = pipeEasySmooth == i
        morphodims = morphoEasy.shape
        #Adjust pipeline cell to morphoeasy size by resizing, then filter so  it's binary again:
        pipecellshort = transform.resize(pipecell[:,:,:], (morphodims), preserve_range= True) 
        overlap = morphoEasy*pipecellshort
        newinold = [] #list with the amount of overlap of the old seed in the new seeds
        candidates = np.unique(overlap)
        for r in candidates:
            if r == 0:
                newinold.append(0) #We skip coincidences with the background
                pass
            else:
                occs = overlap == r
                over = occs.sum()
                newinold.append(over)
        for j in range(len(newinold)):
            if newinold[j] == max(newinold):
                morphonum = candidates[j]
        morphocell = morphoEasy == morphonum
        #Check that tracking is working:
        if i%10 == 0:
            fig, axes = plt.subplots(1,2)
            axes[0].imshow(pipecell[40])
            axes[1].imshow(morphocell[12])
            plt.show
        
        #Great! Now on to measuring volumes!
        
        #We need them to be integer!
        pipecell = pipecell.astype('uint8')
        morphocell = morphocell.astype('uint8')
        Props = measure.regionprops(pipecell)
        pipevols.append(Props[0].area*xyres*xyres*xyres)
        Props2 = measure.regionprops(morphocell)
        morphovols.append(Props2[0].area*xyres*xyres*zres)
        
        #Cell heights:
            
        #FOR PIPELINE
        shell = ndi.binary_dilation(pipecell)
        shell = shell ^ pipecell
        #Now get positions of 1 in shell: They come in the same order as array z, x, y
        coords = np.where(shell==1)
        #Need to check all the values in cells, if value == 0, we add coords to the list and val
        zerocoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if cells3d[coords[0][j], coords[1][j],coords[2][j]] == 0:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                zerocoords = np.vstack((zerocoords, currentcoord))    
        #Thresholding and quality check for values
        thres = filters.threshold_otsu(zerocoords[:,0])
        plt.figure()
        plt.hist(zerocoords[:,0], bins = 20)
        plt.vlines(thres,0, 2000, color = 'red')
        plt.show()
        #Now we can use the threshold to divide the population, put them into two lists
        uppercoords = np.zeros((0,3))
        lowercoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if coords[0][j] > thres:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                uppercoords = np.vstack((uppercoords, currentcoord))
            else:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                lowercoords = np.vstack((lowercoords, currentcoord))
        #make average points from top and bottom and get distance:        
        uppermode = stats.mode(uppercoords[:,0])
        upperpoint = [uppermode[0][0], np.mean(uppercoords[:,1]), np.mean(uppercoords[:,2])]
        lowermode = stats.mode(lowercoords[:,0])
        lowerpoint = [lowermode[0][0], np.mean(lowercoords[:,1]), np.mean(lowercoords[:,2])]
        #Here I could add some quality control, in case the mode is not clearly bigger than the rest of the values...
        #Now distance between two points:
        celllength = abs(uppermode[0][0] - lowermode[0][0])*xyres
        pipelengths.append(celllength)
        
        #FOR MORPHOGRAPH
        shell = ndi.binary_dilation(morphocell)
        shell = shell ^ morphocell
        #Now get positions of 1 in shell: They come in the same order as array z, x, y
        coords = np.where(shell==1)
        #Need to check all the values in cells, if value == 0, we add coords to the list and val
        zerocoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if morphoEasySmooth[coords[0][j], coords[1][j],coords[2][j]] == 0:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                zerocoords = np.vstack((zerocoords, currentcoord))    
        #Thresholding and quality check for values
        thres = filters.threshold_otsu(zerocoords[:,0])
        plt.figure()
        plt.hist(zerocoords[:,0], bins = 20)
        plt.vlines(thres,0, 2000, color = 'red')
        plt.show()
        #Now we can use the threshold to divide the population, put them into two lists
        uppercoords = np.zeros((0,3))
        lowercoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if coords[0][j] > thres:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                uppercoords = np.vstack((uppercoords, currentcoord))
            else:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                lowercoords = np.vstack((lowercoords, currentcoord))
        #make average points from top and bottom and get distance:        
        uppermode = stats.mode(uppercoords[:,0])
        upperpoint = [uppermode[0][0], np.mean(uppercoords[:,1]), np.mean(uppercoords[:,2])]
        lowermode = stats.mode(lowercoords[:,0])
        lowerpoint = [lowermode[0][0], np.mean(lowercoords[:,1]), np.mean(lowercoords[:,2])]
        #Here I could add some quality control, in case the mode is not clearly bigger than the rest of the values...
        #Now distance between two points:
        celllength = abs(uppermode[0][0] - lowermode[0][0])*zres
        morpholengths.append(celllength)
        
            
        #Apical surface
        apicalarea = pipeEasyApical == i
        apicalarea = apicalarea.astype('uint8')
        Props = measure.regionprops(apicalarea)
        if len(Props) == 0:
            pipeapicals.append(0)
        else:
            pipeapicals.append(Props[0].area*xyres*xyres)
        
        morphoapicalarea = morphoEasyApical == morphonum
        morphoapicalarea = morphoapicalarea.astype('uint8')
        Props = measure.regionprops(morphoapicalarea)
        if len(Props) == 0:
            morphoapicals.append(0)
        else:
            morphoapicals.append(Props[0].area*xyres*xyres)
        
        #Basal surface
        basalarea = pipeEasyBasal == i
        basalarea = basalarea.astype('uint8')
        Props = measure.regionprops(basalarea)
        if len(Props) == 0:
            pipebasals.append(0)
        else:
            pipebasals.append(Props[0].area*xyres*xyres)
        
        morphobasalarea = morphoEasyBasal == morphonum
        morphobasalarea = morphobasalarea.astype('uint8')
        Props = measure.regionprops(morphobasalarea)
        if len(Props) == 0:
            morphobasals.append(0)
        else:
            morphobasals.append(Props[0].area*xyres*xyres)


plt.boxplot([pipevols, morphovols])
plt.boxplot([pipeapicals, morphoapicals])
plt.boxplot([pipebasals, morphobasals])
plt.boxplot([pipelengths, morpholengths])


#Now repeat the same for hard segmentation:

        
z, x, y = pipeHardSmooth.shape
pipeHardApical = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if pipeHardSmooth[z-1-c,a,b] == 0:
                pass
            else:
                pipeHardApical[a,b] = pipeHardSmooth[z-1-c,a,b]
                break    
plt.imshow(pipeHardApical)

#And basal:
pipeHardBasal = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if pipeHardSmooth[c,a,b] == 0:
                pass
            else:
                pipeHardBasal[a,b] = pipeHardSmooth[c,a,b]
                break    
plt.imshow(pipeHardBasal) 

#Maybe we have to polish a bit morpho segmentation by eliminating thin slices that were unassigned
morphoHardSmooth = np.zeros_like(morphoHard)
for i in np.unique(morphoHard):
    if i == 0:
        continue
    currentcell = morphoHard == i
    currentcell = currentcell.astype('uint8')
    Props = measure.regionprops(currentcell)
    vol = Props[0].area*xyres*xyres*zres
    if vol>400:
        morphoHardSmooth = np.where(currentcell==1, i, morphoHardSmooth)

#And now apicals and basals:
z, x, y = morphoHardSmooth.shape
morphoHardApical = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if morphoHardSmooth[z-1-c,a,b] == 0:
                pass
            else:
                morphoHardApical[a,b] = morphoHardSmooth[z-1-c,a,b]
                break    
plt.imshow(morphoHardApical)

#And basal:
morphoHardBasal = np.zeros((x,y))
for a in range(y):
    for b in range(x):
        for c in range(z):
            if morphoHardSmooth[c,a,b] == 0:
                pass
            else:
                morphoHardBasal[a,b] = morphoHardSmooth[c,a,b]
                break    
plt.imshow(morphoHardBasal) 

pipevols2 = []
pipeapicals2 = []
pipebasals2 = []
pipelengths2 = []
morphovols2 = []
morphoapicals2 = []
morphobasals2 = []
morpholengths2 = []
for i in np.unique(pipeHardClean):
    if i ==0:
        continue
    else:
        #First, find parallel cell:
        print("on cell " + str(i))
        pipecell = pipeHardSmooth == i
        morphodims = morphoHard.shape
        #Adjust pipeline cell to morphoHard size by resizing, then filter so  it's binary again:
        pipecellshort = transform.resize(pipecell[:,:,:], (morphodims), preserve_range= True) 
        overlap = morphoHard*pipecellshort
        newinold = [] #list with the amount of overlap of the old seed in the new seeds
        candidates = np.unique(overlap)
        for r in candidates:
            if r == 0:
                newinold.append(0) #We skip coincidences with the background
                pass
            else:
                occs = overlap == r
                over = occs.sum()
                newinold.append(over)
        for j in range(len(newinold)):
            if newinold[j] == max(newinold):
                morphonum = candidates[j]
        morphocell = morphoHard == morphonum
        #Check that tracking is working:
        if i%10 == 0:
            fig, axes = plt.subplots(1,2)
            axes[0].imshow(pipecell[40])
            axes[1].imshow(morphocell[12])
            plt.show
        
        #Great! Now on to measuring volumes!
        
        #We need them to be integer!
        pipecell = pipecell.astype('uint8')
        morphocell = morphocell.astype('uint8')
        Props = measure.regionprops(pipecell)
        pipevols2.append(Props[0].area*xyres*xyres*xyres)
        Props2 = measure.regionprops(morphocell)
        morphovols2.append(Props2[0].area*xyres*xyres*zres)
        
        #Cell heights:
            
        #FOR PIPELINE
        shell = ndi.binary_dilation(pipecell)
        shell = shell ^ pipecell
        #Now get positions of 1 in shell: They come in the same order as array z, x, y
        coords = np.where(shell==1)
        #Need to check all the values in cells, if value == 0, we add coords to the list and val
        zerocoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if cells3d[coords[0][j], coords[1][j],coords[2][j]] == 0:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                zerocoords = np.vstack((zerocoords, currentcoord))    
        #Thresholding and quality check for values
        thres = filters.threshold_otsu(zerocoords[:,0])
        plt.figure()
        plt.hist(zerocoords[:,0], bins = 20)
        plt.vlines(thres,0, 2000, color = 'red')
        plt.show()
        #Now we can use the threshold to divide the population, put them into two lists
        uppercoords = np.zeros((0,3))
        lowercoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if coords[0][j] > thres:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                uppercoords = np.vstack((uppercoords, currentcoord))
            else:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                lowercoords = np.vstack((lowercoords, currentcoord))
        #make average points from top and bottom and get distance:        
        uppermode = stats.mode(uppercoords[:,0])
        upperpoint = [uppermode[0][0], np.mean(uppercoords[:,1]), np.mean(uppercoords[:,2])]
        lowermode = stats.mode(lowercoords[:,0])
        lowerpoint = [lowermode[0][0], np.mean(lowercoords[:,1]), np.mean(lowercoords[:,2])]
        #Here I could add some quality control, in case the mode is not clearly bigger than the rest of the values...
        #Now distance between two points:
        celllength = abs(uppermode[0][0] - lowermode[0][0])*xyres
        pipelengths2.append(celllength)
        
        #FOR MORPHOGRAPH
        shell = ndi.binary_dilation(morphocell)
        shell = shell ^ morphocell
        #Now get positions of 1 in shell: They come in the same order as array z, x, y
        coords = np.where(shell==1)
        #Need to check all the values in cells, if value == 0, we add coords to the list and val
        zerocoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if morphoHardSmooth[coords[0][j], coords[1][j],coords[2][j]] == 0:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                zerocoords = np.vstack((zerocoords, currentcoord))    
        #Thresholding and quality check for values
        thres = filters.threshold_otsu(zerocoords[:,0])
        plt.figure()
        plt.hist(zerocoords[:,0], bins = 20)
        plt.vlines(thres,0, 2000, color = 'red')
        plt.show()
        #Now we can use the threshold to divide the population, put them into two lists
        uppercoords = np.zeros((0,3))
        lowercoords = np.zeros((0,3))
        for j in range(len(coords[0])):
            if coords[0][j] > thres:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                uppercoords = np.vstack((uppercoords, currentcoord))
            else:
                currentcoord = np.array((coords[0][j], coords[1][j],coords[2][j]))
                lowercoords = np.vstack((lowercoords, currentcoord))
        #make average points from top and bottom and get distance:        
        uppermode = stats.mode(uppercoords[:,0])
        upperpoint = [uppermode[0][0], np.mean(uppercoords[:,1]), np.mean(uppercoords[:,2])]
        lowermode = stats.mode(lowercoords[:,0])
        lowerpoint = [lowermode[0][0], np.mean(lowercoords[:,1]), np.mean(lowercoords[:,2])]
        #Here I could add some quality control, in case the mode is not clearly bigger than the rest of the values...
        #Now distance between two points:
        celllength = abs(uppermode[0][0] - lowermode[0][0])*zres
        morpholengths2.append(celllength)
        
            
        #Apical surface
        apicalarea = pipeHardApical == i
        apicalarea = apicalarea.astype('uint8')
        Props = measure.regionprops(apicalarea)
        if len(Props) == 0:
            pipeapicals2.append(0)
        else:
            pipeapicals2.append(Props[0].area*xyres*xyres)
        
        morphoapicalarea = morphoHardApical == morphonum
        morphoapicalarea = morphoapicalarea.astype('uint8')
        Props = measure.regionprops(morphoapicalarea)
        if len(Props) == 0:
            morphoapicals2.append(0)
        else:
            morphoapicals2.append(Props[0].area*xyres*xyres)
        
        #Basal surface
        basalarea = pipeHardBasal == i
        basalarea = basalarea.astype('uint8')
        Props = measure.regionprops(basalarea)
        if len(Props) == 0:
            pipebasals2.append(0)
        else:
            pipebasals2.append(Props[0].area*xyres*xyres)
        
        morphobasalarea = morphoHardBasal == morphonum
        morphobasalarea = morphobasalarea.astype('uint8')
        Props = measure.regionprops(morphobasalarea)
        if len(Props) == 0:
            morphobasals2.append(0)
        else:
            morphobasals2.append(Props[0].area*xyres*xyres)


plt.boxplot([pipevols2, morphovols2])
plt.boxplot([pipeapicals2, morphoapicals2])
plt.boxplot([pipebasals2, morphobasals2])
plt.boxplot([pipelengths2, morpholengths2])


#Finall with all:
plt.figure(figsize = (9,8))
plt.subplot(2,2,1)
plt.boxplot([pipevols + pipevols2, morphovols + morphovols2])
plt.xticks([1,2],["Our pipeline", "Morphograph X"])
plt.ylabel("volume (fl)")
plt.subplot(2,2,2)
plt.boxplot([pipeapicals + pipeapicals2, morphoapicals + morphoapicals2])
plt.xticks([1,2],["Our pipeline", "Morphograph X"])
plt.ylabel("apical area (micrometer^2)")
plt.subplot(2,2,3)
plt.boxplot([pipebasals + pipebasals2, morphobasals + morphobasals2])
plt.xticks([1,2],["Our pipeline", "Morphograph X"])
plt.ylabel("basal area (micrometer^2)")
plt.subplot(2,2,4)
plt.boxplot([pipelengths + pipelengths2, morpholengths + morpholengths2])
plt.xticks([1,2],["Our pipeline", "Morphograph X"])
plt.ylabel("cell height (micrometer)")
plt.show()

#Number of cells validated:
allvols = pipevols + pipevols2
allvolsmorph = morphovols + morphovols2

len(allvols)
len(allvolsmorph)
