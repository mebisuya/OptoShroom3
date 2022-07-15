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
                    nslices = nz//5
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
                    bkfactor = bkminval + bstsize*valuelocation[1] -10 #Manual adjustment because with these settings bkground was too stringents and couldnt find a way to fix it
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
                if lap == 2:
                    plt.savefig(graphpath[:-4] + r"SegmentParameters.png", bbox_inches='tight')
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
dirpath = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\20minfiles'
dirpath2 = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hfiles'
dirpath3 = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\c450vfiles'
from os import listdir
filenames = listdir(dirpath)
filepaths = []
for i in range(len(filenames)):
    filepaths.append(join(dirpath, filenames[i]))
filenames2 = listdir(dirpath2)
for i in range(len(filenames2)):
    filepaths.append(join(dirpath2, filenames2[i]))
filenames3 = listdir(dirpath3)
for i in range(len(filenames3)):
    filepaths.append(join(dirpath3, filenames3[i]))

filenamesall = filenames + filenames2 + filenames3
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
    
    if "2h" in filenamesall[q]:
        graphpath = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hgraphs', originalname)
        videopath = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hvideos', originalname)
        datapath =  join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hdata', originalname)
        segpath =  join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hsegments', originalname)
        datapath2 = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\2hdata3class', originalname)
        stimstart = 30/60
        stimend = 150/60
    elif "450" in filenamesall[q]:
        graphpath = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\c450vgraphs', originalname)
        videopath = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\c450vvideos', originalname)
        datapath =  join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\c450vdata', originalname)
        segpath =  join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\c450vsegments', originalname)
        datapath2 = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\c450vdata3class', originalname)
        stimstart = 30/60
        stimend = 150/60
    else:
        graphpath = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\20mingraphs', originalname)
        videopath = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\20minvideos', originalname)
        datapath =  join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\20mindata', originalname)
        segpath =  join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\20minsegments', originalname)
        datapath2 = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\20mindata3class', originalname)
        stimstart = 30/60
        stimend = 50/60

    finalfile, BKfactor, ENfactor = Segment(image[:,:,:,:], sizefilt =minseed3d, lowsize2d =mincell2d, highsize2d = maxcell2d, bigval = biggestval)

    #Reordering the dimensions: 
    
    finalfile2 = np.swapaxes(finalfile, 3,4) # time,z,x,y,chan by time,z,x,chan,y
    finalfile2 = np.swapaxes(finalfile2, 2,3) #time,z,chan,x,y
    imwrite(segpath  + r"3dseg.tif", finalfile2, imagej = True)
    print("seeds saved")
    
    
    #And then get some values on number of cells, number of cells properly segmented, average volume
    #Number of cells
    cellnumber.append(len(np.unique(finalfile2[0,:,3,:,:])))
    kept = 0
    for a in np.unique(finalfile2[0,:,3,:,:]):
        if a in np.unique(finalfile2[-1,:,3,:,:]):
            kept = kept+1
    cellkept.append(kept/len(np.unique(finalfile2[0,:,3,:,:])))
    
    Props = measure.regionprops(finalfile2[0,:,3,:,:])
    vols = []
    for a in range(len(np.unique(finalfile2[0,:,3,:,:]))-1):
        vols.append(Props[a].area*xyres**3)
    cellvolumes.append(np.mean(vols))
    enfactors.append(ENfactor)
    bkfactors.append(BKfactor)
    print("number of cells segmented:" + str(cellnumber))
    print("number of cells tracked:" + str(cellkept))
    print("average volume of nuclei segmented:" + str(cellvolumes))
    
    segments = np.copy(finalfile2[:,:,1,:,:])
    finalfile3 = np.copy(finalfile2[:,:,0:2,:,:])
    frames, slices, chan, xsize, ysize = finalfile2.shape
    for t in range(frames):
        print("smoothing!" + str(t))
        smoothseg = np.zeros_like(finalfile2[t,:,1,:,:])
        cells3d = np.copy(finalfile2[t,:,1,:,:])
        for i in np.unique(finalfile2[:,:,1,:,:]):
            currentcell = cells3d == i
            currentcell = ndi.binary_erosion(currentcell, iterations = 2)
            currentcell = ndi.binary_dilation(currentcell, iterations = 2)
            smoothseg = np.where(currentcell == 1 , i, smoothseg)
        segments[t,:,:,:] = smoothseg
        print("Cleaning...")
        borders = np.zeros_like(segments[0,:,:,:,]) #(t,z,channel1,x,y) 
        borders[0,:,:] = 1
        borders[:,0,:] = 1
        borders[:,:,0] = 1
        borders[-1,:,:] = 1
        borders[:,-1,:] = 1
        borders[:,:,-1] = 1
        #plt.imshow(borders[0:20,20,0:20], interpolation = "none")
        hits = borders*segments[t,:,:,:] #Will take cells in contact with the border!
        listtodelete = np.unique(hits)
        cleanseg = np.zeros_like(segments[t,:,:,:])
        for i in np.unique(segments[t,:,:,:]):
            if i in listtodelete:
                continue
            else:
                cleanseg = np.where(segments[t,:,:,:] == i, i, cleanseg)
        finalfile3[t,:,1,:,:] = cleanseg
    imwrite(segpath + r"3dsmoothandclean.tif", finalfile3, imagej = True)
    
    #Now operate and analyze segmentation
    cells3d = np.copy(finalfile3[:,:,1,:,:])

    #Im going to pick just the cells that are there all the way through
    propercells = np.zeros_like(cells3d)
    for i in np.unique(cells3d[0,:,:,:]):
        if i ==0: 
            continue
        else:
            inall = "yes"
            for t in range(frames):
                if i in cells3d[t,:,:,:]:
                    pass
                else:
                    inall = "no"
            if inall == "yes":  
                print(i)
                propercells = np.where(cells3d ==i, i, propercells)

    #Measuring volumes of cells throught time

    listofcells = []
    listofvols = []
    for i in np.unique(propercells[0,:,:,:]):
        if i == 0:
            continue
        print(i)
        tempvols = []
        for t in range(frames):
            temparray = propercells[t,:,:,:] == i
            temparray = temparray.astype(int)
            if np.amax(temparray) == 0:
                tempvols.append(0)
            else:
                Props = measure.regionprops(temparray)
                tempvols.append(Props[0].area*xyres*xyres*xyres)
        listofvols.append(tempvols)
        listofcells.append(i)
    
    
    #Now on to cell heigth
    #Because apicobasal axis on our monolayers in glass are always aligned through the z axis
    #We can divide the points of a cell in contact with the background in 2 subgroups, apical and basal. Find an average point
    lengths = []
    unilengths = [] #Versions of length counting only z dimension, will always be lower than lengths
    for n in np.unique(propercells[0,:,:,:]):
        if n == 0:
            continue
        print(n)
        templens = []
        tempuni = []
        for t in range(frames):
            temparray = cells3d[t,:,:,:] == n
            shell = ndi.binary_dilation(temparray)
            shell = shell ^ temparray
            #Now get positions of 1 in shell: They come in the same order as array z, x, y
            coords = np.where(shell==1)
            #Need to check all the values in cells, if value == 0, we add coords to the list and val
            zerocoords = np.zeros((0,3))
            for i in range(len(coords[0])):
                if cells3d[t,coords[0][i], coords[1][i],coords[2][i]] == 0:
                    currentcoord = np.array((coords[0][i], coords[1][i],coords[2][i]))
                    zerocoords = np.vstack((zerocoords, currentcoord))    
            #Thresholding and quality check for values
            thres = filters.threshold_otsu(zerocoords[:,0])
            if (t+1)%frames ==0:
                plt.figure()
                plt.hist(zerocoords[:,0], bins = 20)
                plt.vlines(thres,0, 2000, color = 'red')
                plt.show()
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
            templens.append(celllength)
            tempuni.append(abs(uppermode[0][0] - lowermode[0][0]))
        lengths.append(templens)
        unilengths.append(tempuni)
    #Fast idea, for each pixel go down the z, and keep the first value we encounter in the segmentation
    frames, z, y, x = propercells.shape
    apicalcuts = np.zeros((frames,3,x,y))

    for t in range(frames):
        print(t)
        for a in range(y):
            for b in range(x):
                for c in range(z):
                    if finalfile3[t,z-1-c,1,a,b] == 0:
                        pass
                    else:
                        apicalcuts[t,1,a,b] = propercells[t,z-1-c,a,b]
                        apicalcuts[t,2,a,b] = c
                        apicalcuts[t,0,a,b] = finalfile3[t,z-1-c,0,a,b]
                        break

    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(apicalcuts[t,0,:,:])
    plt.subplot(1,3,2)
    plt.imshow(apicalcuts[t,1,:,:])
    plt.subplot(1,3,3)
    plt.imshow(apicalcuts[t,2,:,:])
    plt.show()

    #Save this valuable info as one tif!
    imwrite(segpath + r"apicalcuts.tif", apicalcuts)

    #Let's do the same for basal surfaces!
    basalcuts = np.zeros((frames,3,x,y))

    for t in range(frames):
        print(t)
        for a in range(y):
            for b in range(x):
                for c in range(z):
                    if finalfile3[t,c,1,a,b] == 0:
                        pass
                    else:
                        basalcuts[t,1,a,b] = propercells[t,c,a,b]
                        basalcuts[t,2,a,b] = c
                        basalcuts[t,0,a,b] = finalfile3[t,c,0,a,b]
                        break

    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(basalcuts[t,0,:,:])
    plt.subplot(1,3,2)
    plt.imshow(basalcuts[t,1,:,:])
    plt.subplot(1,3,3)
    plt.imshow(basalcuts[t,2,:,:])
    plt.show()

    #Save this valuable info as one tif!
    imwrite(segpath + r"basalcuts.tif", basalcuts)

    #That worked great!!!! Now we measure apical areas
    areas = []
    for n in np.unique(apicalcuts[0,1,:,:]):
        if n == 0:
            continue
        area = []
        for t in range(frames):
            if n not in apicalcuts[t,1,:,:]:
                #For very weird exceptions when cell does not have apical (likely missegmented)
                area.append(0)
            else:
                tempcell = apicalcuts[t,1,:,:] == n
                tempcell = tempcell.astype(int)
                Props = measure.regionprops(tempcell)
                area.append(Props[0].area*xyres*xyres)
        areas.append(area)

    #and basal areas
    basalareas = []
    for n in np.unique(basalcuts[0,1,:,:]):
        if n == 0:
            continue
        area = []
        for t in range(frames):
            if n not in basalcuts[t,1,:,:]:
                #For very weird exceptions when cell does not have basal (likely missegmented)
                area.append(0)
            else:
                tempcell = basalcuts[t,1,:,:] == n
                tempcell = tempcell.astype(int)
                Props = measure.regionprops(tempcell)
                area.append(Props[0].area*xyres*xyres)
        basalareas.append(area)


    #add more quality filters by eliminating cells that change volume very suddenly
    print(listofcells)

    clearedcells = []
    clearedvols = []
    clearedlengths = []
    clearedunilengths = []
    clearedapicals = []
    clearedbasals = []
    for n in range(len(listofvols)):
        difs = []
        for t in range(len(listofvols[n])):
            if t == 0:
                continue
            else:
                dif = (listofvols[n][t]-listofvols[n][t-1])/listofvols[n][t-1]
                #if the variation in volume is bigger than 50%, the cell is either dividing or mis-segmented
                difs.append(math.sqrt(dif**2))
        if np.amax(difs) > 0.3:
            continue
        else:
            clearedcells.append(listofcells[n])
            clearedvols.append(listofvols[n])
            clearedlengths.append(lengths[n])
            clearedunilengths.append(unilengths[n])
            clearedapicals.append(areas[n])
            clearedbasals.append(basalareas[n])
    print(clearedcells)

    #And now classify cells in and out of stimulation area.  cells in the stimulation area at the beginning of the imaging:
    stimulatedcellspre = np.unique(propercells[0,:,192:320,192:320])
    #This may be suboptimal, I shoud elliminate first cells uncleared

    nonstimulatedcells = []
    stimvols = []
    nonstimvols = []
    stimlengths = []
    stimunilengths = []
    nonstimlengths = []
    nonstimunilengths = []
    stimapicals = []
    nonstimapicals = []
    stimbasals = []
    nonstimbasals = []
    stimulatedcells = []
    for i in range(len(clearedcells)):
        if clearedcells[i] in stimulatedcellspre:
            stimulatedcells.append(clearedcells[i])
            stimvols.append(clearedvols[i])
            stimlengths.append(clearedlengths[i])
            stimunilengths.append(clearedunilengths[i])
            stimapicals.append(clearedapicals[i])
            stimbasals.append(clearedbasals[i])
        else:
            nonstimvols.append(clearedvols[i])
            nonstimlengths.append(clearedlengths[i])
            nonstimunilengths.append(clearedunilengths[i])
            nonstimapicals.append(clearedapicals[i])
            nonstimbasals.append(clearedbasals[i])  
            nonstimulatedcells.append(clearedcells[i])

    #Now all the averages:

    avstimvols = []
    avstimapicals = []
    avstimlengths  = []
    avstimbasals = []
    for t in range(frames):
        temparea = []
        templength = []
        tempvol = []
        tempbasal = []
        for n in range(len(stimulatedcells)):
            temparea.append(stimapicals[n][t])
            templength.append(stimlengths[n][t])
            tempvol.append(stimvols[n][t])
            tempbasal.append(stimbasals[n][t])
        avstimvols.append(np.mean(tempvol))
        avstimapicals.append(np.mean(temparea))
        avstimbasals.append(np.mean(tempbasal))
        avstimlengths.append(np.mean(templength))

    #Nonstimulated cells:
    avnonstimvols = []
    avnonstimapicals = []
    avnonstimbasals = []
    avnonstimlengths  = []
    for t in range(frames):
        temparea = []
        templength = []
        tempvol = []
        tempbasal = []
        for n in range(len(nonstimulatedcells)):
            temparea.append(nonstimapicals[n][t])
            templength.append(nonstimlengths[n][t])
            tempvol.append(nonstimvols[n][t])
            tempbasal.append(nonstimbasals[n][t])
        avnonstimvols.append(np.mean(tempvol))
        avnonstimapicals.append(np.mean(temparea))
        avnonstimlengths.append(np.mean(templength))
        avnonstimbasals.append(np.mean(tempbasal))

    #And now on to the graphs!
    plt.figure(figsize = (9,5))
    plt.title("Evolution of volume. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], stimvols[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avstimvols, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], nonstimvols[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avnonstimvols, c ='magenta', label = "Non Stimulated")
    plt.legend()
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    plt.xlabel("Time (h)")
    plt.ylabel("Volume (femtoliters)")
    #plt.show()
    plt.savefig(graphpath[:-4] + r"CellVolume.png", bbox_inches='tight')

    plt.figure(figsize = (9,5))
    plt.title("Evolution of cell heigth. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], stimlengths[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avstimlengths, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], nonstimlengths[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avnonstimlengths, c ='magenta', label = "Non Stimulated")
    plt.legend()
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    plt.xlabel("Time (h)")
    plt.ylabel("Length in apico-basal axis (um)")
    #plt.show()
    plt.savefig(graphpath[:-4] + r"Cellheigth.png", bbox_inches='tight')

    plt.figure(figsize = (9,5))
    plt.title("Evolution of cell apical area. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], stimapicals[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avstimapicals, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], nonstimapicals[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avnonstimapicals, c ='magenta', label = "Non Stimulated")#0.5 is the time resolution
    plt.legend()
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    plt.xlabel("Time (h)")
    plt.ylabel("Apical area (um2)")
    #plt.show()
    plt.savefig(graphpath[:-4] + r"ApicalArea.png", bbox_inches='tight')

    plt.figure(figsize = (9,5))
    plt.title("Evolution of cell basal area. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], stimbasals[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avstimbasals, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], nonstimbasals[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], avnonstimbasals, c ='magenta', label = "Non Stimulated")#0.5 is the time resolution
    plt.legend()
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    plt.xlabel("Time (h)")
    plt.ylabel("Basal area (um2)")
    #plt.show()
    plt.savefig(graphpath[:-4] + r"BasalArea.png", bbox_inches='tight')

    #Finally I would like to calculate relative valuegraphs: 
    #Volumes: 
    relstimvols = []
    relstimlengths = []
    relstimunilengths = []
    relstimapicals = []
    relstimbasals = []
    for i in range(len(stimvols)):
        tempv = []
        templ = []
        tempunil = []
        tempa = []
        tempb = []
        for t in range(frames):
            tempv.append(stimvols[i][t]/stimvols[i][0])
            templ.append(stimlengths[i][t]/stimlengths[i][0])
            tempunil.append(stimunilengths[i][t]/stimunilengths[i][0])
            tempa.append(stimapicals[i][t]/stimapicals[i][0])
            tempb.append(stimbasals[i][t]/stimbasals[i][0])
        relstimvols.append(tempv)
        relstimlengths.append(templ)
        relstimapicals.append(tempa)
        relstimbasals.append(tempb)
        relstimunilengths.append(tempunil)
        
    relnonstimvols = []
    relnonstimlengths = []
    relnonstimunilengths = []
    relnonstimapicals = []
    relnonstimbasals = []
    for i in range(len(nonstimvols)):
        tempv = []
        templ = []
        tempunil = []
        tempa = []
        tempb = []
        for t in range(frames):
            tempv.append(nonstimvols[i][t]/nonstimvols[i][0])
            templ.append(nonstimlengths[i][t]/nonstimlengths[i][0])
            tempunil.append(nonstimunilengths[i][t]/nonstimunilengths[i][0])
            tempa.append(nonstimapicals[i][t]/nonstimapicals[i][0])
            tempb.append(nonstimbasals[i][t]/nonstimbasals[i][0])
        relnonstimvols.append(tempv)
        relnonstimlengths.append(templ)
        relnonstimunilengths.append(tempunil)
        relnonstimapicals.append(tempa)
        relnonstimbasals.append(tempb)

    #Averages
    relavstimvols = []
    relavstimapicals = []
    relavstimlengths  = []
    relavstimbasals = []
    for t in range(frames):
        temparea = []
        templength = []
        tempvol = []
        tempbasal = []
        for n in range(len(stimulatedcells)):
            temparea.append(relstimapicals[n][t])
            templength.append(relstimlengths[n][t])
            tempvol.append(relstimvols[n][t])
            tempbasal.append(relstimbasals[n][t])
        relavstimvols.append(np.mean(tempvol))
        relavstimapicals.append(np.mean(temparea))
        relavstimbasals.append(np.mean(tempbasal))
        relavstimlengths.append(np.mean(templength))

    #Nonstimulated cells:
    relavnonstimvols = []
    relavnonstimapicals = []
    relavnonstimlengths  = []
    relavnonstimbasals = []
    for t in range(frames):
        temparea = []
        templength = []
        tempvol = []
        tempbasal
        for n in range(len(nonstimulatedcells)):
            temparea.append(relnonstimapicals[n][t])
            templength.append(relnonstimlengths[n][t])
            tempvol.append(relnonstimvols[n][t])
            tempbasal.append(relnonstimbasals[n][t])
        relavnonstimvols.append(np.mean(tempvol))
        relavnonstimapicals.append(np.mean(temparea))
        relavnonstimlengths.append(np.mean(templength))
        relavnonstimbasals.append(np.mean(tempbasal))

    plt.figure(figsize = (9,5))
    plt.title("Evolution of volume. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relstimvols[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavstimvols, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relnonstimvols[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavnonstimvols, c ='magenta', label = "Non Stimulated")
    plt.legend()
    plt.xlabel("Time (h)")
    plt.ylabel("Volume (relative)")
    #plt.show()
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    plt.savefig(graphpath[:-4] + r"CellVolumeRelative.png", bbox_inches='tight')

    plt.figure(figsize = (9,5))
    plt.title("Evolution of cell heigth. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relstimlengths[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavstimlengths, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relnonstimlengths[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavnonstimlengths, c ='magenta', label = "Non Stimulated")
    plt.legend()
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    plt.xlabel("Time (h)")
    plt.ylabel("Length in apico-basal axis (relative)")
    #plt.show()
    plt.savefig(graphpath[:-4] + r"CellheigthRelative.png", bbox_inches='tight')

    plt.figure(figsize = (9,5))
    plt.title("Evolution of cell apical area. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relstimapicals[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavstimapicals, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relnonstimapicals[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavnonstimapicals, c ='magenta', label = "Non Stimulated")#0.5 is the time resolution
    plt.legend()
    plt.xlabel("Time (h)")
    plt.ylabel("Apical area (relative)")
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    #plt.show()
    plt.savefig(graphpath[:-4] + r"ApicalAreaRelative.png", bbox_inches='tight')

    plt.figure(figsize = (9,5))
    plt.title("Evolution of cell basal area. MDCK cells")
    for n in range(len(stimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relstimbasals[n], alpha = 0.2, c ='blue') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavstimbasals, c ='blue', label = "Stimulated")
    for n in range(len(nonstimulatedcells)):
        plt.plot([i*timeres for i in range(frames)], relnonstimbasals[n], alpha = 0.2, c ='magenta') #0.5 is the time resolution
    plt.plot([i*timeres for i in range(frames)], relavnonstimbasals, c ='magenta', label = "Non Stimulated")#0.5 is the time resolution
    plt.legend()
    plt.xlabel("Time (h)")
    plt.ylabel("Basal area (relative)")
    plt.axvspan(stimstart, stimend, color="blue", alpha = 0.1)
    #plt.show()
    plt.savefig(graphpath + r"BasalAreaRelative.png", bbox_inches='tight')


    #Time to export the data! 
    Outputlist = []
    Outputlist.append(stimvols)
    Outputlist.append(nonstimvols)
    Outputlist.append(relstimvols)
    Outputlist.append(relnonstimvols)
    Outputlist.append(stimlengths)
    Outputlist.append(nonstimlengths)
    Outputlist.append(relstimlengths)
    Outputlist.append(relnonstimlengths)
    Outputlist.append(stimapicals)
    Outputlist.append(nonstimapicals)
    Outputlist.append(relstimapicals)
    Outputlist.append(relnonstimapicals)
    Outputlist.append(stimbasals)
    Outputlist.append(nonstimbasals)
    Outputlist.append(relstimbasals)
    Outputlist.append(relnonstimbasals)
    Outputlist.append(stimunilengths)
    Outputlist.append(nonstimunilengths)
    Outputlist.append(relstimunilengths)
    Outputlist.append(relnonstimunilengths)
    DataOut = datapath + "analysis.npy"
    np.save(DataOut, Outputlist)
    print("Data has been saved to" + DataOut)

    #This will be the way to load this data
    #temp = np.load(DataOut, allow_pickle = True)

    #Alternative of analysis by dividing cells between stimulated, adjacent, and non stimulated. Will save everything parallely.

    stimulatedcellspre = np.unique(propercells[4,:,192:320,192:320])
    #To choose adjacent cells, we pick on average the next 2 layers of cells
    #With average apical area of : 140um2, if we assume circle shape:
    # area = pi*r**2 
    #cellradius = math.sqrt(140/math.pi) = 6.7, so diameter  is 13.4, so 2 cells would be 26.8um
    # given that xy resolution is 0.207,
    #that's a total of 134 pixels, I was doing 80, which is a bit more than one layer
    #So maybe I will keep it like that, so the first two rows are selected, but it can always be changed
    adjacentcells = list(np.unique(propercells[4,:,112:192,112:400]))+list(np.unique(propercells[0,:,320:400,112:400]))+list(np.unique(propercells[0,:,192:320,112:192]))+list(np.unique(propercells[0,:,192:320,320:400]))
    #This may be suboptimal, I shoud elliminate first cells uncleared

    listofcells = [[],[],[]] # All going to be triple nested lists,  with #stimulated, #adjacent and #outer cells
    volumes = [[],[],[]] # All going to be triple nested lists,  with #stimulated, #adjacent and #outer cells
    lengths = [[],[],[]] # All going to be triple nested lists,  with #stimulated, #adjacent and #outer cells
    unilengths = [[],[],[]]
    apicals = [[],[],[]] # All going to be triple nested lists,  with #stimulated, #adjacent and #outer cells
    basals = [[],[],[]] # All going to be triple nested lists,  with #stimulated, #adjacent and #outer cells
    for i in range(len(clearedcells)):
        if clearedcells[i] in stimulatedcellspre:
            listofcells[0].append(clearedcells[i])
            volumes[0].append(clearedvols[i])
            lengths[0].append(clearedlengths[i])
            unilengths[0].append(clearedunilengths[i])
            apicals[0].append(clearedapicals[i])
            basals[0].append(clearedbasals[i])
        else:
            if clearedcells[i] in adjacentcells:
                listofcells[1].append(clearedcells[i])
                volumes[1].append(clearedvols[i])
                lengths[1].append(clearedlengths[i])
                unilengths[1].append(clearedunilengths[i])
                apicals[1].append(clearedapicals[i])
                basals[1].append(clearedbasals[i])
            else:
                listofcells[2].append(clearedcells[i])
                volumes[2].append(clearedvols[i])
                lengths[2].append(clearedlengths[i])
                unilengths[2].append(clearedunilengths[i])
                apicals[2].append(clearedapicals[i])
                basals[2].append(clearedbasals[i])
    
    #Sanity check that we have cells in all groups:   
    aretherecells = "yes"
    for i in range(3):
        if len(listofcells[i]) > 0:
            continue
        else:
            aretherecells = "no"
            print("No 3 layer analysis because there are not cells in all 3 classes")
    
    if aretherecells == "yes":
        #Now all the averages:
        def getAverage(nestedlist):
            #For 2 layer nested lists of same lenght == repetitions of same measurements of different experiments
            if len(nestedlist) == 0:
                averagelist = []
            else:
                averagelist = []
                standevup = []
                standevlow = []
                for i in range(len(nestedlist[0])):
                    temp = []
                    for n in range(len(nestedlist)):
                        temp.append(nestedlist[n][i])
                    averagelist.append(np.mean(temp))
            return averagelist
            
    
        avvolumes = [[],[],[]]
        avlengths = [[],[],[]]
        avapicals = [[],[],[]]
        avbasals = [[],[],[]]
        for i in range(3):
            avvolumes[i].append(getAverage(volumes[i]))
            avlengths[i].append(getAverage(lengths[i]))
            avapicals[i].append(getAverage(apicals[i]))
            avbasals[i].append(getAverage(basals[i]))
            
        #And now on to the graphs!
        colorlist = ['blue', 'magenta', 'green']
        labels = ['Stimulated', 'adjacent', 'outercells']
        plt.figure(figsize = (9,5))
        plt.title("Evolution of volume. MDCK cells")
        for a in range(3):
            if len(avvolumes[a]) == 0:
                continue
            for n in range(len(listofcells[a])):
                plt.plot([i*timeres for i in range(frames)], volumes[a][n], alpha = 0.2, c =colorlist[a]) #0.5 is the time resolution
            plt.plot([i*timeres for i in range(frames)], avvolumes[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Volume (femtoliters)")
        #plt.show()
        plt.savefig(graphpath[:-4] + r"CellVolume3Class.png", bbox_inches='tight')
    
        plt.figure(figsize = (9,5))
        plt.title("Evolution of cell heigth. MDCK cells")
        for a in range(3):
            for n in range(len(listofcells[a])):
                plt.plot([i*timeres for i in range(frames)], lengths[a][n], alpha = 0.2, c =colorlist[a]) #0.5 is the time resolution
            plt.plot([i*timeres for i in range(frames)], avlengths[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Length in apico-basal axis (um)")
        #plt.show()
        plt.savefig(graphpath[:-4] + r"Cellheigth3class.png", bbox_inches='tight')
    
        plt.figure(figsize = (9,5))
        plt.title("Evolution of cell apical area. MDCK cells")
        for a in range(3):
            for n in range(len(listofcells[a])):
                plt.plot([i*timeres for i in range(frames)], apicals[a][n], alpha = 0.2, c =colorlist[a]) #0.5 is the time resolution
            plt.plot([i*timeres for i in range(frames)], avapicals[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Apical area (um2)")
        #plt.show()
        plt.savefig(graphpath[:-4] + r"ApicalArea3class.png", bbox_inches='tight')
    
        plt.figure(figsize = (9,5))
        plt.title("Evolution of cell basal area. MDCK cells")
        for a in range(3):
            for n in range(len(listofcells[a])):
                plt.plot([i*timeres for i in range(frames)], basals[a][n], alpha = 0.2, c =colorlist[a]) #0.5 is the time resolution
            plt.plot([i*timeres for i in range(frames)], avbasals[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Basal area (um2)")
        #plt.show()
        plt.savefig(graphpath[:-4] + r"BasalArea3class.png", bbox_inches='tight')
    
        #Finally I would like to calculate relative valuegraphs: 
        #Volumes: 
        relvolumes = [[],[],[]]
        rellengths = [[],[],[]]
        relunilengths = [[],[],[]]
        relapicals = [[],[],[]]
        relbasals = [[],[],[]]
        for n in range(3):
            for i in range(len(listofcells[n])):
                tempv = []
                templ = []
                tempunil = []
                tempa = []
                tempb = []
                for t in range(frames):
                    tempv.append(volumes[n][i][t]/volumes[n][i][0])
                    templ.append(lengths[n][i][t]/lengths[n][i][0])
                    tempunil.append(unilengths[n][i][t]/unilengths[n][i][0])
                    tempa.append(apicals[n][i][t]/apicals[n][i][0])
                    tempb.append(basals[n][i][t]/basals[n][i][0])
                relvolumes[n].append(tempv)
                rellengths[n].append(templ)
                relunilengths[n].append(tempunil)
                relapicals[n].append(tempa)
                relbasals[n].append(tempb)
            
        relavvolumes = [[],[],[]]
        relavlengths = [[],[],[]]
        relavapicals = [[],[],[]]
        relavbasals = [[],[],[]]
        for i in range(3):
            relavvolumes[i].append(getAverage(relvolumes[i]))
            relavlengths[i].append(getAverage(rellengths[i]))
            relavapicals[i].append(getAverage(relapicals[i]))
            relavbasals[i].append(getAverage(relbasals[i]))
    
        plt.figure(figsize = (9,5))
        plt.title("Evolution of volume. MDCK cells")
        for a in range(3):
            for n in range(len(listofcells[a])):
                print(n)
                plt.plot([i*timeres for i in range(frames)], relvolumes[a][n], alpha = 0.2, c =colorlist[a]) 
            plt.plot([i*timeres for i in range(frames)], relavvolumes[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Volume (relative)")
        plt.show()
        plt.savefig(graphpath[:-4] + r"CellVolumeRelative3class.png", bbox_inches='tight')
    
        plt.figure(figsize = (9,5))
        plt.title("Evolution of cell heigth. MDCK cells")
        for a in range(3):
            for n in range(len(listofcells[a])):
                print(n)
                plt.plot([i*timeres for i in range(frames)], rellengths[a][n], alpha = 0.2, c =colorlist[a]) 
            plt.plot([i*timeres for i in range(frames)], relavlengths[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Length in apico-basal axis (relative)")
        #plt.show()
        plt.savefig(graphpath[:-4] + r"CellheigthRelative3class.png", bbox_inches='tight')
    
        plt.figure(figsize = (9,5))
        plt.title("Evolution of cell apical area. MDCK cells")
        for a in range(3):
            for n in range(len(listofcells[a])):
                print(n)
                plt.plot([i*timeres for i in range(frames)], relapicals[a][n], alpha = 0.2, c =colorlist[a]) 
            plt.plot([i*timeres for i in range(frames)], relavapicals[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Apical area (relative)")
        #plt.show()
        plt.savefig(graphpath[:-4] + r"ApicalAreaRelative3class.png", bbox_inches='tight')
    
        plt.figure(figsize = (9,5))
        plt.title("Evolution of cell basal area. MDCK cells")
        for a in range(3):
            for n in range(len(listofcells[a])):
                print(n)
                plt.plot([i*timeres for i in range(frames)], relbasals[a][n], alpha = 0.2, c =colorlist[a]) 
            plt.plot([i*timeres for i in range(frames)], relavbasals[a][0], c =colorlist[a], label = labels[a])
        plt.legend()
        plt.axvspan(stimstart, stimend, color="deepskyblue", alpha = 0.1)
        plt.xlabel("Time (h)")
        plt.ylabel("Basal area (relative)")
        #plt.show()
        plt.savefig(graphpath + r"BasalAreaRelative3class.png", bbox_inches='tight')
    
    
        #Time to export the data! 
        Outputlist = []
        for i in range(3):
            Outputlist.append(volumes[i])
        for i in range(3):
            Outputlist.append(lengths[i])
        for i in range(3):
            Outputlist.append(apicals[i])
        for i in range(3):
            Outputlist.append(basals[i])
        for i in range(3):
            Outputlist.append(relvolumes[i])
        for i in range(3):
            Outputlist.append(rellengths[i])
        for i in range(3):
            Outputlist.append(relapicals[i])
        for i in range(3):
            Outputlist.append(relbasals[i])
        for i in range(3):
            Outputlist.append(unilengths[i])
        for i in range(3):
            Outputlist.append(relunilengths[i])
        DataOut2 = datapath2 + "analysis3class.npy"
        np.save(DataOut2, Outputlist)
        print("Data has been saved to" + DataOut)  

    #We make the 3D array with the cleared cells:
    clearcells = np.zeros_like(propercells)
    for i in clearedcells:
        print(i)
        clearcells = np.where(propercells == i, i, clearcells)

    #Resave the final file, now with the clearcells:
    finalfile3[:,:,1,:,:] = clearcells
    imwrite(segpath+ r"3dsmoothandclean.tif", finalfile3, imagej = True)
    
    
    end = time.time()
    print("time required for image" + str((end - start)/60) + " min")

#Once the screening is done we can plot:
graphpath2 = join(r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\testgraphs', "globalgraph")

plt.figure(figsize = (16,8))
plt.suptitle("")
plt.subplot(1,5,1)
plt.title("Percentage of cells tracked")
plt.boxplot(cellkept)
plt.ylim(0,1.2)
plt.subplot(1,5,2)
plt.title("Number of cells")
plt.boxplot(cellnumber)
plt.ylim(0,100)
plt.subplot(1,5,3)
plt.title("volume of seeds")
plt.boxplot(cellvolumes)
plt.ylim(0,500)
plt.subplot(1,5,4)
plt.title("enfactors")
plt.boxplot(enfactors)
plt.ylim(0,400)
plt.subplot(1,5,5)
plt.title("bkfactors")
plt.boxplot(bkfactors)
plt.ylim(0,400)
plt.savefig(graphpath2 + r"MainValues.png", bbox_inches='tight')
plt.show()

datapath = r'H:\3D segmentation on Virtual Machine\Reversible Stimulation Tests\testdata'

#and save the data!

Outputlist = []
Outputlist.append(cellkept)
Outputlist.append(cellnumber)
Outputlist.append(cellvolumes)
Outputlist.append(enfactors)
Outputlist.append(bkfactors)
DataOut = datapath + "analysis.npy"
np.save(DataOut, Outputlist)
print("Data has been saved to" + DataOut)
    