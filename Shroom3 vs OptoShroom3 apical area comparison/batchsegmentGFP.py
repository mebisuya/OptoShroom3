# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 09:04:49 2022

@author: ara
"""

###Batch 3D segmentation of MDCK stimulated monolayers:

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from skimage import filters, measure, morphology, transform
from skimage.io import imread
from skimage.segmentation import watershed
# importing the image file
from os.path import join
import time
from tifffile import imwrite, TiffFile
import math
from scipy import stats
#Needed subfunctions
def xymask(I, z, enfactor, bkfactor, sigma, topsig):
    masks = np.zeros_like(I)
    enhan = np.zeros_like(I)
    back = np.zeros_like(I)
    for s in range(int(z)):
        #print("We are on slice " + str(s))
        #bkfactor will act as a mask to eliminate background.
        if z < 2:
            imgfin = I
        else:
            imgfin = I[s,:,:] #Assign the stack we are segmenting
        enh = filters.meijering(imgfin, sigmas = np.linspace(0,topsig,4), black_ridges = False)
        enhgauss = ndi.gaussian_filter(enh, sigma)
        np.amax(imgfin)
        enhgauss = enhgauss*np.amax(imgfin)
        inverted = enhgauss < enfactor
        imggaussian = ndi.gaussian_filter(imgfin, sigma)
        #Thresholding 
        foreground = imggaussian>bkfactor
        foreground = ndi.binary_dilation(foreground)
        foreground = ndi.binary_fill_holes(foreground)
        background = np.logical_not(foreground)
        #I don't do erosion on foreground on xy because they are monolayers and there is no background on most  frames
        #foreground = ndi.binary_erosion(foreground, iterations = 5)
        masked = inverted
        masked = ndi.binary_opening(masked)
        masked = ndi.binary_dilation(masked, iterations = 2)
        #masked = ndi.binary_erosion(masked, iterations = 2)
        if z < 2:
            back = background
            enhan = enhgauss
            masks = masked
        else:
            back[s,:,:] = background
            enhan[s,:,:] = enhgauss
            masks[s,:,:] = masked
    return masks, enhan, back

def xzmask(I, y, enfactor, bkfactor, sigma, topsig):
    masks = np.zeros_like(I)
    enhan = np.zeros_like(I)
    back = np.zeros_like(I)
    for s in range(int(y)):
        #print("We are on slice " + str(s))
        #bkfactor will act as a mask to eliminate background.
        if y < 2:
            imgfin = I
        else:
            imgfin = I[:,s,:] #Assign the stack we are segmenting
        enh = filters.meijering(imgfin, sigmas = np.linspace(0,topsig,2), black_ridges = False)
        enhgauss = ndi.gaussian_filter(enh, sigma)
        enhgauss = enhgauss*np.amax(imgfin)
        inverted = enhgauss < enfactor
        imggaussian = ndi.gaussian_filter(imgfin, sigma+1)
        #Thresholding 
        foreground = imggaussian>bkfactor
        #foreground = ndi.binary_dilation(foreground)
        foreground = ndi.binary_fill_holes(foreground)
        foreground = ndi.binary_erosion(foreground, iterations = 5)
        background = np.logical_not(foreground)
        foreground = ndi.binary_erosion(foreground, iterations = 10)
        masked = inverted*foreground
        masked = ndi.binary_dilation(masked, iterations = 2)
        #masked = ndi.binary_erosion(masked, iterations = 2)
        if y < 2:
            back = background
            enhan = enhgauss
            masks = masked
        else:
            back[:,s,:] = background
            enhan[:,s,:] = enhgauss
            masks[:,s,:] = masked
    return masks, enhan, back

def yzmask(I, x, enfactor, bkfactor, sigma, topsig):
    masks = np.zeros_like(I)
    enhan = np.zeros_like(I)
    back = np.zeros_like(I)
    for s in range(int(x)):
        #print("We are on slice " + str(s))
        #bkfactor will act as a mask to eliminate background.
        if x < 2:
            imgfin = I
        else:
            imgfin = I[:,:,s] #Assign the stack we are segmenting
        enh = filters.meijering(imgfin, sigmas = np.linspace(0,topsig,2), black_ridges = False)
        enhgauss = ndi.gaussian_filter(enh, sigma)
        enhgauss = enhgauss*np.amax(imgfin)
        inverted = enhgauss < enfactor
        imggaussian = ndi.gaussian_filter(imgfin, sigma+1) #Increased this sigma because we only want to separate betwee bcground and fground
        #Thresholding 
        foreground = imggaussian>bkfactor
        #foreground = ndi.binary_dilation(foreground)
        foreground = ndi.binary_fill_holes(foreground)
        foreground = ndi.binary_erosion(foreground, iterations = 5)
        background = np.logical_not(foreground)
        foreground = ndi.binary_erosion(foreground, iterations = 10)
        masked = inverted*foreground
        masked = ndi.binary_dilation(masked, iterations = 2)
        #masked = ndi.binary_erosion(masked, iterations = 2)

        if x < 2:
            back = background
            enhan = enhgauss
            masks = masked
        else:
            back[:,:,s] = background
            enhan[:,:,s] = enhgauss
            masks[:,:,s] = masked
    return masks, enhan, back

from numba import njit, prange

@njit(parallel=True)
def finalize(triple, trp, confidence):
    z, x, y = triple.shape
    for a in range(z):
        #print(a)
        for b in range(x):
            #print(b)
            for c in prange(y):
                #print(c)
                vals = np.unique(trp[a,b,c,:]) #Count the values including ALL neighbours
                numbvals = vals.shape
                #print(numbvals[0])
                if numbvals[0] > 1:
                    if a > z-2 or b > x-2 or c > y-2:
                        triple[a,b,c] = trp[a,b,c,0]
                        confidence[a,b,c] = 1
                    else:
                        neib = np.append(trp[a,b,c,:],trp[a+1,b,c,:])
                        neib = np.append(neib, trp[a-1,b,c,:])
                        neib = np.append(neib, trp[a,b+1,c,:])
                        neib = np.append(neib, trp[a,b-1,c,:])
                        neib = np.append(neib, trp[a,b,c+1,:])
                        neib = np.append(neib, trp[a,b,c-1,:])
                        neib = np.append(neib, trp[a-1,b-1,c-1,:])
                        neib = np.append(neib, trp[a-1,b-1,c,:])
                        neib = np.append(neib, trp[a-1,b-1,c+1,:])
                        neib = np.append(neib, trp[a-1,b+1,c-1,:])
                        neib = np.append(neib, trp[a-1,b+1,c,:])
                        neib = np.append(neib, trp[a-1,b+1,c+1,:])
                        neib = np.append(neib, trp[a+1,b+1,c-1,:])
                        neib = np.append(neib, trp[a+1,b+1,c,:])
                        neib = np.append(neib, trp[a+1,b+1,c+1,:])
                        neib = np.append(neib, trp[a-1,b+1,c-1,:])
                        neib = np.append(neib, trp[a-1,b+1,c,:])
                        neib = np.append(neib, trp[a-1,b+1,c+1,:])
                        neib = np.append(neib, trp[a-1,b,c-1,:])
                        neib = np.append(neib, trp[a+1,b,c-1,:])
                        neib = np.append(neib, trp[a-1,b,c+1,:])
                        neib = np.append(neib, trp[a+1,b,c+1,:])
                        neib = np.append(neib, trp[a,b-1,c-1,:])
                        neib = np.append(neib, trp[a,b+1,c-1,:])
                        neib = np.append(neib, trp[a,b-1,c+1,:])
                        neib = np.append(neib, trp[a,b+1,c+1,:])
                        #Can simplify this with a loop!
                        vals = np.unique(neib)
                        counts = np.zeros((0))
                        #print(neib)
                        if numbvals[0] == 2:
                            confidence[a,b,c] = 2
                        else:
                            confidence[a,b,c] = 3
                        for d in vals:
                            occ = neib == d
                            tot = occ.sum()
                            counts = np.append(counts, tot)
                        #I have no choice but to iterate again through the list
                        #print(max(counts))
                        sp = counts.shape
                        for e in range(sp[0]):
                            if counts[e]==max(counts):
                                triple[a,b,c] = vals[e]
                                break
                            else:
                                pass
                else:
                    confidence[a,b,c] = 1
                    counts = np.zeros((0))
                    for d in vals:
                        occ = trp[a,b,c,:] == d
                        tot = occ.sum()
                        counts = np.append(counts, tot)
                    #I have no choice but to iterate again through the list
                    #print(max(counts))
                    sp = counts.shape
                    for e in range(sp[0]):
                        if counts[e]==max(counts):
                            triple[a,b,c] = vals[e]
                            break
                        else:
                            pass

def plotScore2D(array, enlist, bklist, lap):
    #given an [n,n,4] dimensional array, in which the last 4 dimensions are number of cells, area, background factor, and score
    x, y, dims = array.shape
    columns = dims//2
    rows = dims//2
    plt.figure(figsize = (columns*7, rows*7))
    plt.suptitle("Factor selection, iteration " + str(lap))
    nplot = 1
    for n in range(dims): 
        plt.subplot(rows, columns, nplot)
        nplot = nplot+1
        if n == 0: 
            plt.title("cell number")
            plt.imshow(array[:,:,n])
        elif n == 1:
            plt.title("area")
            plt.imshow(array[:,:,n])
        elif n == 2:
            plt.title("area stand dev")
            plt.imshow(array[:,:,n])
        else:
            plt.title("Score")
            plt.imshow(array[:,:,n], cmap = 'plasma')
        plt.xticks(np.arange(len(bklist)), labels = bklist)
        plt.yticks(np.arange(len(enlist)), labels = enlist)
        plt.xlabel("background factor")
        plt.ylabel("Enhanced factor")
    plt.show()
    return


def Segment(I, sizefilt = 250, lowsize2d = 400, highsize2d = 1000, divfact = 1, sigma = 2.5, bkfactor = 0.55, enfactor = 0.55, ensig =4, bigval = 65536):
    #We start by reshaping the volume to homogeneously sized voxels
    if len(I.shape) < 4:
        stack, x, y = I.shape
        tpoints = 1
    else:
        tpoints, stack, x, y = I.shape
    ratio = zres/xyres
    z = int(stack*ratio) #Calculates the number of Z stacks we should have approximately if x-y were the same resolution as z
    #New dimensions for reduced segmentation
    nx = x//divfact
    ny = y//divfact
    nz = z//divfact #In case we want to reduce the stacks for faster processing
    segmented = np.zeros((tpoints,nz,nx,ny,4))
    num = 1 # This is a variable to count cells, every time we add a new valid seed we will count it as new unless we trace it...
    for t in range(tpoints):
        print("We are on timepoint " + str(t))
        if tpoints <2:
            img = transform.resize(I[:,:,:], (nz,nx,ny), preserve_range= True) 
        else:
            img = transform.resize(I[t,:,:,:], (nz,nx,ny), preserve_range= True)    
        if t == 0:
            exslice = int(nz//2)

            print("Looking for optimal parameters for segmentation")
            #Original images to get thresholds, from vertical, which will surely have some background
            #We will use two threshold filters, one to detect the background to not background. 
            resteps = 10 #Number of steps in each direction
            bkmaxval = bigval//10 #maximum values to get tested, minimum always 0
            bkminval = 0
            enmaxval = bigval//10
            enminval = 0 
            nslices = 5
            for lap in range(3): #3 rounds of iteration to find the best value
                bstsize = (bkmaxval-bkminval)/resteps
                estsize = (enmaxval-enminval)/resteps
                screen = np.zeros((resteps,resteps,4))
                bklist = []
                enlist = []
                if lap == 1:
                    nslices = 7
                    #We add more slices to fine tune!
                for b in range(resteps):
                    print(100*b/resteps)
                    bkfactor = bkminval + b*bstsize
                    bklist.append(np.around(bkfactor,3))
                    for e in range(resteps):
                        enfactor = enminval + e*estsize
                        if b == 0:
                            enlist.append(np.around(enfactor, 3))
                        masklist = []
                        for i in range(nslices-1):
                            maskedxz, enhancedxz, gaussianxz= xzmask(img[:,(i+1)*ny//nslices,:], 1, enfactor, bkfactor, sigma, ensig)
                            maskedyz, enhancedyz, gaussianyz= yzmask(img[:,:,(i+1)*nx//nslices], 1, enfactor, bkfactor, sigma, ensig)
                            maskedxy, enhancedxy, gaussianxy= xymask(img[(i+1)*nz//nslices,:,:], 1, enfactor, bkfactor, sigma, ensig)
                            masklist.append(maskedyz)
                            masklist.append(maskedxz)
                            masklist.append(maskedxy)
                        #We will make measuremements on Masked xy
                        scorelist = []
                        arealist =[]
                        sdlist = []
                        ncellslist = []
                        for mask in masklist:
                            labelled, ncellstotal = ndi.label(mask)
                            labelled = labelled.astype(int)
                            Props = measure.regionprops(labelled)
                            if ncellstotal>0:
                                areas = []
                                areasforsd = []
                                for n in range(ncellstotal):
                                    areasforsd.append(Props[n].area)
                                    if Props[n].area < lowsize2d:
                                        continue
                                    else:
                                        if Props[n].area < highsize2d:
                                            areas.append(Props[n].area)
                                        else:
                                            pass
                                if len(areas) == 0:
                                    areas.append(0)
                                area = np.mean(areas)
                                sd = np.std(areasforsd)
                                ncells = len(areas) #Remake ncells after filtering small areas
                            else:
                                area = 0
                                sd = 0
                                ncells=0
                            if sd == 0:
                                score = 0
                            else:
                                score = ncells*ncells/ncellstotal #it's the number of cells that made it through the filter
                            arealist.append(ncellstotal)
                            sdlist.append(sd)
                            ncellslist.append(ncells)
                            scorelist.append(score)
                        screen[e,b,0] = np.mean(ncellslist)
                        screen[e,b,1] = np.mean(arealist)
                        screen[e,b,2] = np.mean(sdlist)
                        screen[e,b,3] = np.mean(scorelist)
                        
                        
                #Gaussian filter to score to reduce outliers
                screen[:,:,3] = ndi.filters.gaussian_filter(screen[:,:,3], sigma = 0.5)
                plotScore2D(screen, enlist, bklist, lap)
                tuplelocation = np.where(screen[:,:,3] == np.amax(screen[:,:,3]))
                valuelocation = [tuplelocation[0][len(tuplelocation[0])//2], tuplelocation[1][len(tuplelocation[1])//2]]
                if lap == 2:
                    bkfactor = bkminval + bstsize*valuelocation[1] +30 #Manual adjustment because with these settings bkground was too soft
                else:
                    bkfactor = bkminval + bstsize*valuelocation[1]
                enfactor = enminval + estsize*valuelocation[0]
                print("Selected bkfactor = " + str(np.around(bkfactor,3)) + ", enfactor = " + str(np.around(enfactor, 3)) + " and enh sigma = " + str(np.around(ensig,3)))
                print("score was = " + str(np.amax(screen[:,:,3])))
                #Now let's show how it looks with the picked factors!
                maskedxy, enhancedxy, backgroundxy= xymask(img[exslice,:,:], 1, enfactor, bkfactor, sigma, ensig)
                maskedyz, enhancedyz, backgroundyz= yzmask(img[:,ny//2,:], 1, enfactor, bkfactor, sigma, ensig)
                
                
                plt.figure()
                plt.suptitle("bkfactor =" + str(np.around(bkfactor,3)) + ", enhfactor = " + str(np.around(enfactor,3)) + " and ensigma = " + str(np.around(ensig,3)))
                plt.subplot(2,3,1)
                plt.title("Enhanced")
                plt.imshow(enhancedxy, interpolation = None)
                plt.subplot(2,3,2)
                plt.title("Background")
                plt.imshow(backgroundxy, interpolation = None)
                plt.subplot(2,3,3)
                plt.title("Seeds")
                plt.imshow(maskedxy, interpolation = None)
                plt.subplot(2,3,4)
                plt.imshow(enhancedyz, interpolation = None)
                plt.subplot(2,3,5)
                plt.imshow(backgroundyz, interpolation = None)
                plt.subplot(2,3,6)
                plt.imshow(maskedyz, interpolation = None)
                plt.show()
                
                
                bkminval = np.around(bkfactor - bstsize*1, 3)
                bkmaxval = np.around(bkfactor + bstsize*1, 3)
                enminval =  np.around(enfactor - estsize*1, 3)
                enmaxval = np.around(enfactor + estsize*1, 3)
        print("Segmenting seeds in the 3 axes")
        mxy, exy, bxy= xymask(img, nz, enfactor, bkfactor, sigma, ensig)
        print("z done")
        mxz, exz, bxz= xzmask(img, ny, enfactor, bkfactor, sigma, ensig)
        print("y done")
        myz, eyz, byz= yzmask(img, nx, enfactor, bkfactor, sigma, ensig)
        print("x done")
        
        print("Polishing seeds")
        
        mult1 = mxy*mxz
        mult2 = mult1*myz        
        
        #same for backgound seed
        back1 = bxz
        back2 = back1*byz #Now we keep it for watershed!
        
        
        for s in range(int(nz)):
            #print(s)
            mult2[s,:,:] = ndi.binary_dilation(mult2[s,:,:])
            mult2[s,:,:] = ndi.binary_fill_holes(mult2[s,:,:])
            mult2[s,:,:] = ndi.binary_erosion(mult2[s,:,:], iterations = 3)
        print("z done") 
        for s in range(int(ny)):
            #print(s)
            mult2[:,s,:] = ndi.binary_dilation(mult2[:,s,:])
            mult2[:,s,:] = ndi.binary_fill_holes(mult2[:,s,:])
            mult2[:,s,:] = ndi.binary_erosion(mult2[:,s,:], iterations = 3)
        print("y done")
        for s in range(int(nx)):
            #print(s)
            mult2[:,:,s] = ndi.binary_dilation(mult2[:,:,s])
            mult2[:,:,s] = ndi.binary_fill_holes(mult2[:,:,s])
            mult2[:,:,s] = ndi.binary_erosion(mult2[:,:,s], iterations = 3)
        print("x done")
        
        
        labelled, ncells = ndi.label(mult2)
        print("Filtering seeds")
        filtered = np.zeros_like(img)
        for i in np.unique(labelled):
            if i == 0:
                pass
            else:
                ocurrences = labelled == i # TRUE only in cells where volume has that value
                vol = ocurrences.sum()
                #print(vol)
                if vol > sizefilt:
                    filtered= np.where(labelled == i, num, filtered)
                    num = num+1

        zred = int(nz)
        yred = int(ny)
        xred = int(nx)
        
        #Here I am including the tracing based on the overlap of the seeds.
        #I need to obtain the overlaps from a different file than the one I'm modifying to avoid errors
        if t > 0:
            filtered2 = np.copy(filtered)
            print("Tracking seeds")
            for i in np.unique(segmented[t-1,:,:,:,3]):
                if i == 0:
                    pass
                else:
                    ocurrences = segmented[t-1,:,:,:,3] == i
                    overlap = ocurrences*filtered2 #Array only with the seeds overlaping
                    newinold = [] #list with the amount of overlap of the old seed in the new seeds
                    candidates = np.unique(overlap)
                    for r in candidates:
                        if r == 0:
                            newinold.append(0) #We skip coincidences with the background
                            pass
                        else:
                            occs = overlap == r
                            over = occs.sum()
                            newinold.append(over) # We collect all the pixels in common with the candidate old seeds
                    if len(newinold) ==1:
                        pass
                    else:
                        for j in range(len(newinold)):
                            if newinold[j] == max(newinold): #For the old seed that shares more pixels, we make sure that it's also it's max overlaping seed
                                filtered = np.where(filtered2 == candidates[j], i, filtered)
                            else:
                                pass
            #A final pass to make all the numbers go one after the other and not leaving gaps:
            oldseeds = np.amax(np.unique(segmented[t-1,:,:,:,3]))
            for i in np.unique(filtered):
                if i > np.amax(segmented[t-1,:,:,:,3]):
                    oldseeds = oldseeds + 1
                    filtered = np.where(filtered == i, oldseeds, filtered)
                else:
                    pass
            num = oldseeds + 1 #preparing num for the next timepoint
        
        print("Watersheding")
        tripleWater = np.zeros((zred,xred,yred,3)) # One for each dimension!
        #On xy
        #We add the background seed!
        filtered = filtered + back2*2000 #2000 is going to be the number for background
        for s in range(zred):
            #print("We are on slice " + str(s))
            water = watershed(exy[s,:,:], filtered[s,:,:])
            tripleWater[s,:,:,0] = water
        print("z done")
        for s in range(xred):
            #print("We are on slice " + str(s))
            water = watershed(exz[:,s,:], filtered[:,s,:])
            tripleWater[:,s,:,1] = water
        print("x done")
        for s in range(yred):
            #print("We are on slice " + str(s))
            water = watershed(eyz[:,:,s], filtered[:,:,s])
            tripleWater[:,:,s,2] = water
        print("y done")
        
        finalsegment = np.zeros_like(tripleWater[:,:,:,0])
        confidence =  np.zeros_like(tripleWater[:,:,:,0])
        print("finalizing")
        finalize(finalsegment, tripleWater, confidence)
        
        #Now we can remove the background value, from now on it will be 0
        finalsegment = np.where(finalsegment > 1999, 0, finalsegment)
        filtered = np.where(filtered > 1999, 0, filtered)
        #And finally we store them in a multichannel stack
        segmented[t,:,:,:,0] = img
        segmented[t,:,:,:,1] = finalsegment
        segmented[t,:,:,:,2] = confidence
        segmented[t,:,:,:,3] = filtered
        print("timepoint done")
    segmented = segmented.astype('int16') #it turns resized image into 0s !
    return segmented, bkfactor, enfactor

#Import the file
 #Taking filename without .tif to use it when saving data
dirpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\Shroom3VolumeMeasurements\GFPCAAXFiles'
from os import listdir
filenames = listdir(dirpath)
filepaths = []
for i in range(len(filenames)):
    filepaths.append(join(dirpath, filenames[i]))
filenamesall = filenames
enfactors = []
bkfactors = []
cellnumber = []
cellkept = []
cellvolumes = []
for q in range(len(filepaths)):
    start = time.time()
    filepath = filepaths[q]

    with TiffFile(filepath) as tif:
        image = tif.asarray()
        imagej_metadata = tif.imagej_metadata
    print('We are on file ' + filenamesall[q])
    originalname = filenamesall[q][:-4]
    #Let's resize the image for proper 3D segmentation, making squared voxels. 
    xyres = 0.2071602 #xy um/pixel I used the same xy res for all the analysis
    zres = float(imagej_metadata['spacing'])  #z um/pixel
    print("zres is " + str(zres))
    timeres = 0.166666666 #10 mins/60mins ( in hours)
    biggestval = np.amax(image)
    #The area pir**2 == 314um2 one slice of a cell. Very small cells should be 40um2 minimum then:
    mincell2d = 10/(xyres**2) #We will use this to pick the filters!
    maxcell2d =  250/(xyres**2) #Cells must not be bigger than this value 
    #Also MDCK cells are on average between 1000 and 2000 um3 in volume.
    #"The spherical nucleus typically occupies about 10 percent of a eukaryotic cell's volume"
    #Then a nuclei should minimum 100um3. If we go to a extreme, let's say 40um3, we will discard all seeds below that
    minseed3d = 20/(xyres**3) #Because we will have already rescaled the image, we use only xyres
    datapath =  join(r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\Shroom3VolumeMeasurements\GFPData', originalname)
    segpath =  join(r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\Shroom3VolumeMeasurements\GFPCAAXSegments', originalname)


    finalfile, BKfactor, ENfactor = Segment(image[:,0,:,:], sizefilt =minseed3d, lowsize2d =mincell2d, highsize2d = maxcell2d, bigval = biggestval)

    #Reordering the dimensions: 
    
    finalfile2 = np.swapaxes(finalfile, 3,4) # z,x,y,chan by time,z,x,chan,y
    finalfile2 = np.swapaxes(finalfile2, 2,3) #z,chan,x,y
    imwrite(segpath  + r"3dseg.tif", finalfile2, imagej = True)
    print("seeds saved")
    
    
    #And then get some values on number of cells, number of cells properly segmented, average volume
    #Number of cells
    
    segments = np.copy(finalfile2[0,:,1,:,:])
    finalfile3 = np.copy(finalfile2[0,:,0:2,:,:])
    frames, slices, chan, xsize, ysize = finalfile2.shape
    print("smoothing!")
    smoothseg = np.zeros_like(finalfile2[0,:,1,:,:])
    cells3d = np.copy(finalfile2[0,:,1,:,:])
    for i in np.unique(finalfile2[:,:,1,:,:]):
        currentcell = cells3d == i
        currentcell = ndi.binary_erosion(currentcell, iterations = 2)
        currentcell = ndi.binary_dilation(currentcell, iterations = 2)
        smoothseg = np.where(currentcell == 1 , i, smoothseg)
    segments[:,:,:] = smoothseg
    print("Cleaning...")
    borders = np.zeros_like(segments) #(t,z,channel1,x,y) 
    borders[0,:,:] = 1
    borders[:,0,:] = 1
    borders[:,:,0] = 1
    borders[-1,:,:] = 1
    borders[:,-1,:] = 1
    borders[:,:,-1] = 1

    hits = borders*segments[:,:,:] #Will take cells in contact with the border!
    listtodelete = np.unique(hits)
    cleanseg = np.zeros_like(segments[:,:,:])
    for i in np.unique(segments[:,:,:]):
        if i in listtodelete:
            continue
        else:
            cleanseg = np.where(segments[:,:,:] == i, i, cleanseg)
    finalfile3[:,1,:,:] = cleanseg
    imwrite(segpath + r"3dsmoothandclean.tif", finalfile3, imagej = True)
    
    #Now operate and analyze segmentation
    cells3d = np.copy(finalfile3[:,1,:,:])
    propercells = cells3d

    #Measuring volumes of cells throught time

    listofcells = []
    listofvols = []
    for i in np.unique(propercells[:,:,:]):
        if i == 0:
            continue
        print(i)
        temparray = propercells == i
        temparray = temparray.astype(int)
        if np.amax(temparray) == 0:
            tempvols.append(0)
        else:
            Props = measure.regionprops(temparray)
            listofvols.append(Props[0].area*xyres*xyres*xyres)
            listofcells.append(i)
    
    print(listofvols)
    #Now on to cell heigth
    #Because apicobasal axis on our monolayers in glass are always aligned through the z axis
    #We can divide the points of a cell in contact with the background in 2 subgroups, apical and basal. Find an average point
    lengths = []
    unilengths = [] #Versions of length counting only z dimension, will always be lower than lengths
    for n in np.unique(propercells):
        if n == 0:
            continue
        print(n)
        temparray = cells3d[:,:,:] == n
        shell = ndi.binary_dilation(temparray)
        shell = shell ^ temparray
        #Now get positions of 1 in shell: They come in the same order as array z, x, y
        coords = np.where(shell==1)
        #Need to check all the values in cells, if value == 0, we add coords to the list and val
        zerocoords = np.zeros((0,3))
        for i in range(len(coords[0])):
            if cells3d[coords[0][i], coords[1][i],coords[2][i]] == 0:
                currentcoord = np.array((coords[0][i], coords[1][i],coords[2][i]))
                zerocoords = np.vstack((zerocoords, currentcoord))    
        #Thresholding and quality check for values
        thres = filters.threshold_otsu(zerocoords[:,0])
        #Now we can use the threshold to divide the population, put them into two lists
        uppercoords = np.zeros((0,3))
        lowercoords = np.zeros((0,3))
        for i in range(len(coords[0])):
            if coords[0][i] > thres:
                currentcoord = np.array((coords[0][i], coords[1][i],coords[2][i]))
                uppercoords = np.vstack((uppercoords, currentcoord))
            else:
                currentcoord = np.array((coords[0][i], coords[1][i],coords[2][i]))
                lowercoords = np.vstack((lowercoords, currentcoord))
        #make average points from top and bottom and get distance:        
        uppermode = stats.mode(uppercoords[:,0])
        upperpoint = [uppermode[0][0], np.mean(uppercoords[:,1]), np.mean(uppercoords[:,2])]
        lowermode = stats.mode(lowercoords[:,0])
        lowerpoint = [lowermode[0][0], np.mean(lowercoords[:,1]), np.mean(lowercoords[:,2])]
        #Here I could add some quality control, in case the mode is not clearly bigger than the rest of the values...
        #Now distance between two points:
        celllength = math.sqrt(((upperpoint[0]-lowerpoint[0])**2)+((upperpoint[1]-lowerpoint[1])**2)+((upperpoint[2]-lowerpoint[2])**2))*xyres

        lengths.append(celllength)
        unilengths.append(abs(uppermode[0][0] - lowermode[0][0])*xyres)
    
    
    #Fast idea, for each pixel go down the z, and keep the first value we encounter in the segmentation
    z, y, x = propercells.shape
    apicalcuts = np.zeros((3,x,y))
    for a in range(y):
        for b in range(x):
            for c in range(z):
                if finalfile3[z-1-c,1,a,b] == 0:
                    pass
                else:
                    apicalcuts[1,a,b] = propercells[z-1-c,a,b]
                    apicalcuts[2,a,b] = c
                    apicalcuts[0,a,b] = finalfile3[z-1-c,0,a,b]
                    break

    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(apicalcuts[0,:,:])
    plt.subplot(1,3,2)
    plt.imshow(apicalcuts[1,:,:])
    plt.subplot(1,3,3)
    plt.imshow(apicalcuts[2,:,:])
    plt.show()

    #Save this valuable info as one tif!
    imwrite(segpath + r"apicalcuts.tif", apicalcuts)

    #Let's do the same for basal surfaces!
    basalcuts = np.zeros((3,x,y))

    for a in range(y):
        for b in range(x):
            for c in range(z):
                if finalfile3[c,1,a,b] == 0:
                    pass
                else:
                    basalcuts[1,a,b] = propercells[c,a,b]
                    basalcuts[2,a,b] = c
                    basalcuts[0,a,b] = finalfile3[c,0,a,b]
                    break

    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(basalcuts[0,:,:])
    plt.subplot(1,3,2)
    plt.imshow(basalcuts[1,:,:])
    plt.subplot(1,3,3)
    plt.imshow(basalcuts[2,:,:])
    plt.show()

    #Save this valuable info as one tif!
    imwrite(segpath + r"basalcuts.tif", basalcuts)

    #That worked great!!!! Now we measure apical areas
    areas = []
    for n in np.unique(apicalcuts[1,:,:]):
        if n == 0:
            continue
        else:
            tempcell = apicalcuts[1,:,:] == n
            tempcell = tempcell.astype(int)
            Props = measure.regionprops(tempcell)
            areas.append(Props[0].area*xyres*xyres)

    #and basal areas
    basalareas = []
    for n in np.unique(basalcuts[1,:,:]):
        if n == 0:
            continue
        else:
            tempcell = basalcuts[1,:,:] == n
            tempcell = tempcell.astype(int)
            Props = measure.regionprops(tempcell)
            basalareas.append(Props[0].area*xyres*xyres)
    chosen = np.unique(apicalcuts[1,192:320,192:320])
    chosenvols = []
    chosenaps = []
    chosenbas = []
    chosenlens = []
    for i in range(len(listofcells)):
        if listofcells[i] in chosen:
            chosenvols.append(listofvols[i])
            chosenaps.append(areas[i])
            chosenbas.append(basalareas[i])
            chosenlens.append(unilengths[i])
    

    
    #Time to export the data! 
    Outputlist = []
    Outputlist.append(chosenvols)
    Outputlist.append(chosenaps)
    Outputlist.append(chosenbas)
    Outputlist.append(chosenlens)
    Outputlist.append(chosen)
    DataOut = datapath + "analysis.npy"
    np.save(DataOut, Outputlist)
    print("Data has been saved to" + DataOut)
    
    end = time.time()
    print("time required for image" + str((end - start)/60) + " min")


    