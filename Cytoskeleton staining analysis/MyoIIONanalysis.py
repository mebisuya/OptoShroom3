# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:43:06 2022

@author: ara
"""

##Code to analyze Myosin stainings

#The plan is to load every image, select apical slice to make a max projectiom
#with 3 slices. Then segment on the gfp channel, and use the mask to get the
#average value per pixel inside and outside of the mask. 

#We will need:

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from skimage import filters, measure, morphology, transform, filters
from skimage.io import imread
from tifffile import imwrite, TiffFile
from os import listdir
from os.path import join

path = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\cytoskeleton analysis\MyoIIb ON'
graphpath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\cytoskeleton analysis\graphs'
datapath = r'\\ebisuya.embl.es\ebisuya\Guillermo\Optoshroom3 paper\Python scripts and raw data\cytoskeleton analysis\MyoIIbONdata'

filenames = listdir(path)

for i in range(len(filenames)):
    filepath = join(path, filenames[i])
    with TiffFile(filepath) as tif:
        image = tif.asarray()
        imagej_metadata = tif.imagej_metadata
    z, chan, y ,x = image.shape #z, chan, y, x #GFP IS ON 1!
    #We first need the user to define what is apical:
    plt.figure(figsize = (20,20))
    for a in range(15):
        plt.subplot(3,5,a+1)
        plt.title(str(z-5-a))
        plt.imshow(image[-5-a,1,:,:])
    plt.show()
    apical = input("Where is apical?")
    apical = int(apical)
    gfp = np.amax(image[apical-1:apical+2,1,:,:], axis = 0)
    myo = np.amax(image[apical-1:apical+2,0,:,:], axis = 0)
    
    
    #Now we make an  mask with top 5% intensity
    flatgfp = gfp.flatten()
    sortgfp = np.sort(flatgfp)
    filt = gfp>(sortgfp[-13107])-1 #top 5%
    

    
    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(gfp)
    plt.subplot(1,3,2)
    plt.imshow(myo)
    plt.subplot(1,3,3)
    plt.imshow(filt)
    plt.savefig(graphpath  + "\\" + filenames[i][:-4]+  r"segment.png", bbox_inches='tight')
    plt.show()
    
    #And now we make the measurement:
    onvals = []
    offvals = []
    for n in range(y):
        for d in range(x):
            if filt[n,d] == True:
                onvals.append(myo[n,d])
            else:
                offvals.append(myo[n,d])
    
    onavg = np.mean(onvals)
    offavg = np.mean(offvals)
    ratio = onavg/offavg
    
    #Now we can export!!
    

    
    Outputlist = []
    Outputlist.append(onavg)
    Outputlist.append(offavg)
    Outputlist.append(ratio)
    Outputlist.append(onvals)
    Outputlist.append(offvals)
    DataOut = datapath + "\\" + filenames[i][:-4]+ "analysis.npy"
    np.save(DataOut, Outputlist)

        

