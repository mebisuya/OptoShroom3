# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 09:38:36 2021

@author: ara
"""
# New version of folding analysis, only with the functions to generate graphs of axis skeleton and foldings.
# We will measure curvatures again to check results are the same and I din't change anything important during measurement process

#Copy pasted from previous version

  
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from skimage import filters, measure, morphology
import math
#import tifffile
#import seaborn as sns

#First let's import image and define path
filename = r'20220327_Sample7.tif'
path = r'C:\Users\ara\Documents\Python Scripts\OptoShroom3 folding matrigel beads\Samples'
DataOutputFolder = r"C:\Users\ara\Documents\Python Scripts\OptoShroom3 folding matrigel beads\foldingoutput"
GraphOutputFolder = r"C:\Users\ara\Documents\Python Scripts\OptoShroom3 folding matrigel beads\Graphs\foldinggraphs"
from os.path import join

#Images only 1 by 1 to avoid saturating the system
filepath = join(path, filename)
I = io.imread(filepath)

print("file has been loaded")
SampleClass = input("Is this sample a control (C) or induced (I)") #Define if sample is Control (C) or Induced 

ImageInUse = np.copy(I[:,:,0,:,:])
frames = ImageInUse.shape[0]

def ProjectAndSegment(x):
    Maxproject = np.max(x[:,:,:,:], axis=1) 
    #This function should work if dimensions of ALL files are [t,z,x,y,channel] and channels in same order!! 
    #Now segmentation, has to be specific to frame number of each file
    frames = Maxproject.shape[0]
    masks = np.zeros_like(Maxproject)
    indexs = np.zeros(masks.shape[0], dtype=int)
    labelled_all = np.zeros_like(Maxproject)
    biggestMask = np.zeros_like(Maxproject)
    for t in range(frames):
        otsu = filters.threshold_otsu(Maxproject[t,:,:])
        mask = Maxproject[t,:,:] > otsu * 0.6
        mask = ndi.binary_opening(mask)
        mask = ndi.binary_closing(mask, iterations=3)
        mask = ndi.binary_fill_holes(mask)
        masks[t,:,:] = mask
    #Now that we have the mask, will have to select biggest one. To do so we label and get biggest area 
        labelled_all[t,:,:], nr = ndi.label(masks[t, :,:])
        props = measure.regionprops(labelled_all[t,:,:])
        area_t = np.array([prop['area'] for prop in props])
        if t == 0: #30thMarch2020 Added this to track and keep always the colony that is biggest at t=0
            biggest = max(enumerate(area_t), key=lambda x : x[1])[0] + 1
            indexs[t] = biggest
            biggestMask[t,:,:] = labelled_all[t,:,:] == biggest
            #plt.figure()
            #plt.imshow(masks[0,:,:])
        else: #When not t=0 we find the colony with biggest match
            temporal = np.zeros([0,2],dtype=int)
            for i in np.unique(labelled_all[t,:,:]):
                if i != 0: #excluding  0 to avoid getting the background!
                    currentcol = labelled_all[t,:,:]==i
                    test1 = np.logical_and(currentcol, biggestMask[t-1,:,:]) #To track to previous colony
                    if True in test1:
                        testint=test1.astype(int)
                        Props=measure.regionprops(testint)
                        x = Props[0].area
                        temporal = np.append(temporal, [[i,x]], axis=0)
                    else:
                        continue
                else:
                    continue
            maxim=np.amax(temporal, axis=(0,1))
            d = np.where(temporal==maxim)
            row=int(d[0])
            realcol=temporal[row, 0]
            biggestMask[t,:,:] = labelled_all[t,:,:] == realcol
    return Maxproject, biggestMask

#Function defined!! 
#Now we can implement to all the colonies and get projection and biggest mask systematically
print('Doing projection')
projections, mainmasks = ProjectAndSegment(ImageInUse)
print('Projection done!')

ProjAreas = []
for t in range(frames):    
    Masksprops = measure.regionprops(mainmasks[t,:,:])
    InitialArea = np.array([prop['area'] for prop in Masksprops])
    InitialArea = InitialArea[0]
    ProjAreas.append(InitialArea)
print("Areas have been measured")

#Ok, objectives are: find main axes rotate main file according to it,
#Run 2D segmentation, reslice, segmentation of reslice. For now not worth it
#To run 3D segmentation 

#Now we have projections and 2D masks! Let's find axes and rotate
#We will need some functions

def RadToDeg(a):
    b=a*180/math.pi
    if b< 0:
        b=360 + b
    else:
        pass
    return b


def getPerimeter(mask):
    Rotperimeter = np.zeros_like(mask)
    Rotperimeter = ndi.binary_dilation(mask) ^ mask
    return Rotperimeter


def FindMainAngle(mainmask):
    arraylength = np.zeros((0,2))
    for i in range(360):
        temprot = ndi.interpolation.rotate(mainmask[0,:,:], i, axes=(0,1))
        labelled, nr = ndi.label(temprot[:,:])
        props = measure.regionprops(labelled)
        Centroid = props[0]['centroid'] #we always want this at timepoint 0
        CentroidY = int(Centroid[0])
        tempperim = getPerimeter(temprot)
        linetocheck = tempperim[CentroidY,:]
        pointposition = np.where(linetocheck == 1)
        pointposition = pointposition[0]
        distance = []
        if pointposition.shape[0] > 0:
            for d in range(pointposition.shape[0]-1):
                c= d+1
                distance.append(pointposition[c]-pointposition[d])
            e=[i, max(distance)]
            arraylength = np.vstack((arraylength, e))
        else:
            pass
    maxvalue=np.amax(arraylength[:,1], axis=0)
    mainangles=arraylength[np.where(arraylength[:,1] == maxvalue)]
    # We can get several angles where distance is the same.
    #We will pick the one in the middle. ?
    pick = mainangles.shape[0]//2
    mainangle = mainangles[pick, 0]
    return mainangle

#Now we find angle and rotate:
print("finding main angle, may take a while...")
angle = FindMainAngle(mainmasks)
print("Got the angle, rotating files")
rotMask = ndi.interpolation.rotate(mainmasks, angle, axes=(1,2))
rotProj = ndi.interpolation.rotate(projections, angle, axes=(1,2))
rotFile = ndi.interpolation.rotate(ImageInUse, angle, axes=(2,3))
print("Files rotated!")


#Now let's get main and minor axis values in rotated files:
def getMainAxis(mainmask): #New version taking main axis through centroid too!
    labelled, nr = ndi.label(mainmask[:,:])
    props = measure.regionprops(labelled)
    Centroid = props[0]['centroid'] #we always want this at timepoint 0
    CentroidY = int(Centroid[0])
    tempperim = getPerimeter(mainmask)
    linetocheck = tempperim[CentroidY,:]
    pointposition = np.where(linetocheck == 1)
    pointposition = pointposition[0]
    distance = []
    if pointposition.shape[0] > 0:
        for d in range(pointposition.shape[0]-1):
            c= d+1
            distance.append(pointposition[c]-pointposition[d])
        MainAxis = CentroidY 
        maxvalue = max(distance)
    else:
        pass
    return MainAxis, maxvalue

print("Getiing main axis")
MainAxis, MainSize = getMainAxis(rotMask[0,:,:])

#This could be improved to either get the perimeter once or get it every time
def getMinorAxis(Perimeter, mainaxis):
    #We put perimeter into coordinate points
    pointposition = np.where(Perimeter == 1)
    positions=np.c_[pointposition[0], pointposition[1]]
    r = mainaxis
        #Compute their distance of the points on each row, and keep the maximum of it.
    line = positions[np.where(positions[:,0] == r)]
    if line.shape[0] > 1:
        distance = (line[-1,1]-line[0,1]) #First and last point!
    else:
        print("error, perimeter in main axis has 1 point or less")
    mid = distance/2
    minoraxis = line[0,1] + mid
    minoraxis=int(minoraxis) 
    # This is the minor axis, now let's get it's size
    differences=np.zeros((0, 2))
    lineY = positions[np.where(positions[:,1] == minoraxis)]
    distanceminor = []
    if lineY.shape[0] > 1:
        for d in range(lineY.shape[0]-1):
            c= d+1
            distanceminor.append(lineY[c,0]-lineY[d,0])
        e=[minoraxis, max(distanceminor)]
        differences=np.vstack((differences, e))
    else:
        print("error, perimeter in minor axis has 1 point or less")
    #Here we have the biggest differences for every column. 
    maxvalue=np.amax(differences[:,1], axis=0)
    return minoraxis, maxvalue

print("Getiing minor axis")
perimeter0 = getPerimeter(rotMask[0,:,:])
MinorAxis, MinorSize = getMinorAxis(perimeter0, MainAxis)

AspectRatio = MainSize/MinorSize
print("Aspect ratio of this colony is " + str(AspectRatio))
#Control Check
#print("Your colony and axes look like this")


def getAxisSizes(Mask2D, MainAxis, minoraxis):
    frames = Mask2D.shape[0]
    MainSizes = []
    MinorSizes = []
    for t in range(frames):
        mask = Mask2D[t,:,:]
        Perimeter = np.zeros_like(mask)
        Perimeter = ndi.binary_dilation(mask) ^ mask
    #We put perimeter into coordinate points
        pointposition = np.where(Perimeter == 1)
        positions=np.c_[pointposition[0], pointposition[1]]
        distance = []
        line = positions[np.where(positions[:,0] == MainAxis)]
        if line.shape[0] > 1:
            for d in range(line.shape[0]-1):
                c= d+1
                distance.append(line[c,1]-line[d,1])
                Maxtempmain = max(distance)
            else:
                pass
        MainSizes.append(Maxtempmain)
        #NowMinorAxis:
        lineY = positions[np.where(positions[:,1] == minoraxis)]
        distanceminor = []
        if lineY.shape[0] > 1:
            for d in range(lineY.shape[0]-1):
                c= d+1
                distanceminor.append(lineY[c,0]-lineY[d,0])
            maxtempminor = max(distanceminor)
        else:
            print("error, perimeter in minor axis has 1 point or less")
        MinorSizes.append(maxtempminor)
    return MainSizes, MinorSizes

print('getting all axes sizes')
MainAxSizes, MinorAxSizes = getAxisSizes(rotMask, MainAxis, MinorAxis)

#Cool here we already have some starting data. Now on to lateral curvature


def ResliceMainAxis2(file, mainaxis):
    top = mainaxis + 3
    down = mainaxis - 3
    frames = file.shape[0]
    stacks = file.shape[1]
    top_reslice = np.zeros([frames, stacks, \
                        np.shape(file)[3]], dtype='uint16') #[frames, z, x]
    for t in range(frames):
        top_reslice[t, :,:] = np.max(file[t,::-1, down:top,:], axis=1)
        #We are inverting z to flip upside down
    return top_reslice        

def ResliceMinorAxis2(file, minoraxis):
    top = minoraxis + 3
    down = minoraxis - 3
    frames = file.shape[0]
    stacks = file.shape[1]
    top_reslice = np.zeros([frames, stacks, \
                        np.shape(file)[2]], dtype='uint16') #[frames, z, y]
    for t in range(frames):
        top_reslice[t,:,:] = np.max(file[t,::-1, :, down:top], axis=2)
    return top_reslice

print('Reslicing...')
#test = ResliceMainAxis(Masks3D[4], MainAxisList[4])
#plt.imshow(test[7,:,:], interpolation='none')
MainAxesReslicefile = ResliceMainAxis2(rotFile, MainAxis)
MinorAxesReslicefile = ResliceMinorAxis2(rotFile, MinorAxis)

#plt.imshow(MainAxesReslicefile[12,:,:], interpolation='none')
#plt.imshow(MinorAxesReslicefile[12,:,:], interpolation='none')

#Now I need to crop the areas with Pixel 0. *Appeared with rotation. 
#This way I-ll avoid them from compromising segmentation.
#I'll project on the Y axis, and get rid of all planes in which the value is 0
def CropTheFrame(MainMax, MinorMax):
    frames = MainMax.shape[0]
    yaxis = MainMax.shape[1] #This X and Y are subjective! Not from the original image
    yaxis2 = MinorMax.shape[1] 
    xaxis = MainMax.shape[2] 
    xaxis2 = MinorMax.shape[2] 
    ProjYMain = np.max(MainMax[0,:,:], axis=0) #Makin a line with Max value on y.
    ProjYMinor = np.max(MinorMax[0,:,:], axis=0)
    ColsAre0Main = np.where(ProjYMain == 0)
    ColsAre0Minor =np.where(ProjYMinor == 0)
    ColsAre0Main = ColsAre0Main[0]
    ColsAre0Minor = ColsAre0Minor[0]
    FinalColsMain = xaxis - len(ColsAre0Main)
    FinalColsMinor = xaxis2 - len(ColsAre0Minor)
    CroppedMain = np.zeros(shape=(frames, yaxis, FinalColsMain))
    CroppedMinor = np.zeros(shape=(frames, yaxis2, FinalColsMinor))    
    a = 0
    b = 0 
    for x in range(xaxis):
        if ProjYMain[x] > 0:
            CroppedMain[:,:, a] = MainMax[:,:,x]
            a = a + 1
        else:
            pass
    for x in range(xaxis2):
        if ProjYMinor[x] > 0:
            CroppedMinor[:,:, b] = MinorMax[:,:,x]
            b = b + 1
        else:
            pass
    return CroppedMain, CroppedMinor



Maincrop, Minorcrop = CropTheFrame (MainAxesReslicefile, MinorAxesReslicefile)

#Great! Now we are ready to segment, get skeleton and measure curvature. 

def SegmentAProjection(Proj):
    frames = Proj.shape[0]
    masks = np.zeros_like(Proj)
    indexs = np.zeros(masks.shape[0], dtype=int)
    labelled_all = np.zeros_like(Proj)
    biggestMask = np.zeros_like(Proj)
    for t in range(frames):
        gaussian = ndi.filters.gaussian_filter(Proj[t,:,:], sigma = 1)
        triangle = filters.threshold_triangle(gaussian)
        mask = gaussian > triangle * 1
        mask = ndi.binary_dilation(mask)
        mask = ndi.binary_closing(mask, iterations=3)
        mask = ndi.binary_opening(mask)
        mask = ndi.binary_fill_holes(mask)
        masks[t,:,:] = mask
    #Now that we have the mask, will have to select biggest one. To do so we label and get biggest area 
        labelled_all[t,:,:], nr = ndi.label(masks[t,:,:])
        labelled_all = labelled_all.astype(int)
        props = measure.regionprops(labelled_all[t,:,:])
        area_t = np.array([prop['area'] for prop in props])
        if t == 0: #30thMarch2020 Added this to track and keep always the colony that is biggest at t=0
            biggest = max(enumerate(area_t), key=lambda x : x[1])[0] + 1
            indexs[t] = biggest
            biggestMask[t,:,:] = labelled_all[t,:,:] == biggest
        else: #When not t=0 we find the colony with biggest match
            temporal = np.zeros([1,2],dtype=int)
            for i in np.unique(labelled_all[t,:,:]):
                if i != 0: #excluding  0 to avoid getting the background!
                    currentcol = labelled_all[t,:,:]==i
                    test1 = np.logical_and(currentcol, biggestMask[0,:,:])
                    if True in test1:
                        testint=test1.astype(int)
                        Props=measure.regionprops(testint)
                        x = Props[0].area
                        temporal = np.append(temporal, [[i,x]], axis=0)
                    else:
                        continue
                else:
                    continue
            maxim=np.amax(temporal, axis=(0,1))
            d = np.where(temporal==maxim)
            if d[0].shape[0] > 1:
                row = int(d[0][0])
            else:
                row=int(d[0])
            realcol=temporal[row, 0]
            biggestMask[t,:,:] = labelled_all[t,:,:] == realcol
    return biggestMask

print('Segmenting')
resmainmask = SegmentAProjection(Maincrop)
resminormask = SegmentAProjection(Minorcrop)

#plt.imshow(resminormask[12,:,:], interpolation='none')
#plt.imshow(Minorcrop[12,:,:], interpolation='none')


#Nice! let's skeletonize and get data!

def findSkeletonv2(reslice, sigma=10):
    sks = np.zeros_like(reslice)
    frames = np.shape(reslice)[0]
    for t in range(frames):
        if np.max(reslice[t,:,:]) > 0:
                gauss = filters.gaussian(reslice[t,:,:], sigma)
                th = filters.threshold_otsu(gauss)
                sks[t,:,:] = morphology.skeletonize(gauss > th)
    return(sks)

print('Getting skeletons...')
mainsks = findSkeletonv2(resmainmask)
minorsks = findSkeletonv2(resminormask)


##Cool time to measure!

def PointsInOrder(skeletons):
    frames = np.shape(skeletons)[0]
    IndexesMain = np.where(skeletons == 1)
    #With np.where we get two arrays, values per row (y) that satisfy the condition and then per column. So whe have y,x
    IndexMain =np.column_stack((IndexesMain[0], IndexesMain[1],IndexesMain[2]))
    orderedpoints=np.zeros((0, 5))
    for t in range(frames):
        print(t)
        a=IndexMain[np.where(IndexMain[:,0] == t)]
        #Picking point indexes on current timepoint
        if np.shape(a)[0] < 3:
            print(r"DANGER!No skeleton found on timepoint" + str(t))
            continue
        else:
            pass
        
        nbnumber=np.zeros((0, 4)) #Will be [time, z, x, number of neighbours] index not needed because is kept on the first for iteration
        #First a loop to check number of neighbours, this will help us find the first one and add the next neighbour. 
        for i in range(a.shape[0]):
            x=a[i,2]
            z=a[i,1]
            #Here we test every possible neighbour and put them in a list
            #Neighbours 
            upleft = np.where((a[:,2] == x-1) & (a[:,1] == z+1) )
            up = np.where((a[:,2] == x) & (a[:,1] == z+1) )
            upright = np.where((a[:,2] == x+1) & (a[:,1] == z+1) )
            left = np.where((a[:,2] == x-1) & (a[:,1] == z) )
            right = np.where((a[:,2] == x+1) & (a[:,1] == z) )
            downleft = np.where((a[:,2] == x-1) & (a[:,1] == z-1) )
            down =np.where((a[:,2] == x) & (a[:,1] == z-1)) 
            downright = np.where((a[:,2] == x+1) & (a[:,1] == z-1))
            #We put them all together in one array, number of rows is number of neighbours, and we have the row index for every neighbour. 
            neighbourscheck=np.concatenate((upleft[0], up[0], upright[0], left[0], right[0], downleft[0], down[0], downright[0]), axis=None)
            thisround=np.array(([t, z, x, neighbourscheck.size]))
            nbnumber=np.vstack((nbnumber, thisround))
        #Now we have the list with number of neighbours, we select the first most to the left (neighbour equals 1)
        b=nbnumber[np.where(nbnumber[:,3] == 1)]
        if b.shape[0] > 2:
            print("Branched skeleton! Will skip smallest branch") # Added 6th april 2020
            #Eliminate smallest skeleton Added 4th February 2021
            minorx = np.amin(b[:,2])
            firstpoint = b[np.where(b[:,2] == minorx),:] # We pick the as starting point the one of the selected that is most to the left (smaller x)
            firstpoint = firstpoint[0]
            firstpoint = firstpoint[0]
        elif  b.shape[0] < 2:
            print("Danger, only one single-neighbour point!")
            #This is probably because one of the ends is curled (has 2 or 3 neighbours)
            #That means a point nearby must have 3 neighbours
            b = np.vstack((b, nbnumber[np.where(nbnumber[:,3] == 3)])) # We add them so we can start from that one!
            minorx = np.amin(b[:,2])
            firstpoint = b[np.where(b[:,2] == minorx),:] # We pick the as starting point the one of the selected that is most to the left (smaller x)
            firstpoint = firstpoint[0]
            firstpoint = firstpoint[0] # Eliminating unnecesary array layers
        else: #When points with 1 neighbour are only 2
            if b[0,2]>b[1,2]:
                firstpoint=b[1,:]
            else:
                firstpoint=b[0,:]
        #We have the starting point, lets find it in our nbnumber to start putting them ordered in the bigfile
        c=np.where((nbnumber[:,2] == firstpoint[2]) & (nbnumber[:,1] == firstpoint[1]))
        c=np.concatenate((c))
        c=c.item() #Little trick to just get a scalar value and not an array which is always more messy
        n=0
        old = [] #This is going to be the list of already added points, so we don't go back or repeat!
        for d in range(a.shape[0]):
            # c is the position of the actual point (starting by left limit). We add the  point to the list, we find the next one and discard the older one. 
            # n is the number in the order, we start with 0 
            orderedpoints=np.vstack((orderedpoints, [nbnumber[c,0] , n, nbnumber[c,1], nbnumber[c,2], nbnumber[c,3]]))
            n=n+1
            #Using code from the other loop, find neighbours: In this case, doesn't matter if we use a or nbnumber, both have time,z and x in the same order
            x=nbnumber[c,2]
            z=nbnumber[c,1]
            upleft = np.where((a[:,2] == x-1) & (a[:,1] == z+1) )
            up = np.where((a[:,2] == x) & (a[:,1] == z+1) )
            upright = np.where((a[:,2] == x+1) & (a[:,1] == z+1) )
            left = np.where((a[:,2] == x-1) & (a[:,1] == z) )
            right = np.where((a[:,2] == x+1) & (a[:,1] == z) )
            downleft = np.where((a[:,2] == x-1) & (a[:,1] == z-1) )
            down =np.where((a[:,2] == x) & (a[:,1] == z-1)) 
            downright = np.where((a[:,2] == x+1) & (a[:,1] == z-1))
            #We put them all together in one array, number of rows is number of neighbours, and we have the row index for every neighbour. 
            neighbourscheck=np.concatenate((upleft[0], up[0], upright[0], left[0], right[0], downleft[0], down[0], downright[0]), axis=None)
            # In almost all cases here we will get positions of 2 neigbours, to make sure we pick the next one and not one that is already on the list
            #But it must also work if there is only one neighbour
            if np.size(neighbourscheck) == 1:
                if len(old) != 0:
                    print("we are on a dead end")
                    break;
                else:
                    old.append(c)
                    c=neighbourscheck.item()
            elif np.size(neighbourscheck) > 2:
                #We are at a branching point, evaluating branch lengths...
                old.append(c)
                candidates = [] #List of candidate paths!
                for r in range(len(neighbourscheck)):
                    if neighbourscheck[r] in old:
                        pass
                    else:
                        candidates.append(neighbourscheck[r])
                #Now we need to find if c1 or c2 leads to the bigger branch
                if len(candidates) == 0:
                    print("All neighbours are old, we reached to an end")
                    break;
                else:
                    candidatelength = [] #List to add the length of each possible branch!
                    for branch in range(len(candidates)):
                        c1t = candidates[branch]#Temporal point to iterate through
                        oldt = old.copy()
                        oldt.append(c)
                        c1length = 0
                        for l in range(a.shape[0]):
                            x=nbnumber[c1t,2]
                            z=nbnumber[c1t,1]
                            upleft = np.where((a[:,2] == x-1) & (a[:,1] == z+1) )
                            up = np.where((a[:,2] == x) & (a[:,1] == z+1) )
                            upright = np.where((a[:,2] == x+1) & (a[:,1] == z+1) )
                            left = np.where((a[:,2] == x-1) & (a[:,1] == z) )
                            right = np.where((a[:,2] == x+1) & (a[:,1] == z) )
                            downleft = np.where((a[:,2] == x-1) & (a[:,1] == z-1) )
                            down =np.where((a[:,2] == x) & (a[:,1] == z-1)) 
                            downright = np.where((a[:,2] == x+1) & (a[:,1] == z-1))
                            neighbourscheck=np.concatenate((upleft[0], up[0], upright[0], left[0], right[0], downleft[0], down[0], downright[0]), axis=None)
                            if np.size(neighbourscheck) == 1:
                                break
                            else:
                                for r in range(len(neighbourscheck)):
                                    if neighbourscheck[r] in old:
                                        pass
                                    else:
                                        oldt.append(c1t) #To simpliify we add the first one that is not old
                                        c1t = neighbourscheck[r]
                                        break;
                            c1length += 1
                        candidatelength.append(c1length)
                #print("tail 2 is " + str(c2length))
                #position of longest tail
                #print(candidatelength)
                pmax = np.where(candidatelength == np.amax(candidatelength))
                pmax = pmax[0]
                c = candidates[pmax[0]]# We choose the value that gave to the longest tale
 
            else:
                if neighbourscheck[0] in old:
                    old.append(c)
                    c = neighbourscheck[1]
                else:
                    old.append(c)
                    c = neighbourscheck[0]
        print("This skeleton is of size:")
        print(len(orderedpoints[np.where(orderedpoints[:,0] == t)]))
    return orderedpoints

print('Ordering points, getting curvature')
mainpoints = PointsInOrder(mainsks)
minorpoints = PointsInOrder(minorsks)

##Next, go on to measurements!

#fast test
  
    
def FindCurvature(points, Csize, Csep):
    Circlesize = Csize #Distance between the points of the circle
    Circlesep = Csep #How often we pick points for circles, with 1 we pick every pixel for a circle with its following 2 points. 
    frames = np.amax(points[:,0], axis=0) +1
    Curvatures = np.zeros((0, 3))
    for t in range(int(frames)):
        Pointstemp = points[np.where(points[:,0] == t)]
        CurvaturesTemp=np.zeros((0, 3)) #Will store [orderpoint of the first point of the circle (cs), curvature ]
        for i in range(int((Pointstemp.shape[0]-2*Circlesize)/Circlesep)):  #This way it will stop at the first point of the last circle
            cs=i*Circlesep  #Circle starts at
            pointA=np.array([Pointstemp[cs,2], Pointstemp[cs,3]])
            pointB=np.array([Pointstemp[cs+Circlesize,2], Pointstemp[cs+Circlesize,3]])
            pointC=np.array([Pointstemp[cs+Circlesize*2,2], Pointstemp[cs+Circlesize*2,3]])        
            sidea=math.sqrt(((pointA[0]-pointB[0])**2)+((pointA[1]-pointB[1])**2))
            sideb=math.sqrt(((pointB[0]-pointC[0])**2)+((pointB[1]-pointC[1])**2))
            sidec=math.sqrt(((pointC[0]-pointA[0])**2)+((pointC[1]-pointA[1])**2))
            s=(sidea+sideb+sidec)/2
            if s*(s-sidea)*(s-sideb)*(s-sidec) <= 0: #This happens when the three dots are in a line, which makes curvature == 0
                Curv = 0
            else:
                K=math.sqrt(s*(s-sidea)*(s-sideb)*(s-sidec))
                Curv=4*K/(sidea*sideb*sidec)
            #Adding positive and negative 3rd April2020
            AB = [pointA[0]-pointB[0], pointA[1]-pointB[1]]
            BC = [pointB[0]-pointC[0], pointB[1]-pointC[1]]
            angleABC = math.atan2(BC[1], BC[0]) - math.atan2(AB[1], AB[0])
            angleABC = RadToDeg(angleABC)
            if angleABC >= 180: # if right turn, negative curvature
                Curv = -Curv
            else:
                pass
            CurvaturesTemp=np.vstack((CurvaturesTemp, [t, cs, Curv]))
        Curvatures = np.vstack((Curvatures, CurvaturesTemp))
    return Curvatures

#CIRCLE SIZE MUST BE WAY SMALLER THAN NUMBER OF POINTS IN SMALLER SKELETON

currentnum = []
for t in range(int(frames)):
    tempMain = mainpoints[np.where(mainpoints[:,0] == t)]
    tempMinor = minorpoints[np.where(minorpoints[:,0] == t)]
    ShpMain = np.shape(tempMain)[0]
    ShpMinor = np.shape(tempMinor)[0]
    currentnum.append(ShpMain)
    currentnum.append(ShpMinor)

print('smallest size of skeleton is: ')
print(np.amin(currentnum))

print('Current circle size is 20')

maincurv = FindCurvature(mainpoints, 20, 1) #I'm choosing step size ==20
minorcurv = FindCurvature(minorpoints, 20, 1)
print("curvatures have been measured")
#Until here, we would export this together with ax size data. To be analyzed all in group

#But let's see how average curvature behaved:

time = np.linspace(0, frames-1, frames)
time = time.astype(int)
realtime = time*2

mainmean = []
minormean = []
for i in range(frames):
    Mainstemp = maincurv[np.where(maincurv[:,0] == i)]
    Minorstemp = minorcurv[np.where(minorcurv[:,0] == i)]
    MainMeant = np.mean(Mainstemp[:,2], axis=0)
    MinorMeant = np.mean(Minorstemp[:,2], axis=0)
    mainmean.append(MainMeant)
    minormean.append(MinorMeant)

print("mean curvatures have been calculated")


#Now let's save skeletons, axes and measurements in graphs
from matplotlib.animation import FuncAnimation, PillowWriter  
from matplotlib.colors import LinearSegmentedColormap

#Custom color map for green fluorescence
cdict = {'red':   ((0.0,  0.0, 0.0),
                   (1.0,  0.0, 0.0)),

         'green': ((0.0,  0.0, 0.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (1.0,  0.0, 0.0))}

fluogreen = LinearSegmentedColormap('fluogreen', cdict)

print('we are saving projections')
fig = plt.figure(figsize=(10,10))
plt.axis('off')
fig.tight_layout(pad=0)
def update(tp):  
    MainImage =np.zeros_like(rotMask[tp,:,:])
    MainImage[MainAxis-1:MainAxis+1,:] = 1
    MinorImage = np.zeros_like(rotMask[tp,:,:])
    MinorImage[:,MinorAxis-1:MinorAxis+1] = 1
    MainMask = np.ma.masked_where(MainImage<1, MainImage)
    MinorMask = np.ma.masked_where(MinorImage<1, MinorImage)
    plt.imshow(rotProj[tp,:,:], interpolation='none', cmap= fluogreen)
    #plt.imshow(rotMask[tp,:,:], interpolation='none', alpha=0.2) #We exclude mask for now
    plt.imshow(MainMask, interpolation='none', cmap='Greys')
    plt.imshow(MinorMask, interpolation='none', cmap='Greys')
    
ani = FuncAnimation(fig, update, frames = time)  

writer = PillowWriter(fps=10)  
ani.save(GraphOutputFolder + r"\\" + filename[:-4]+ "projectionsm2.gif", writer=writer) 


print('we are saving main axis of image')
fig = plt.figure(figsize=(10,2))
plt.axis('off')
fig.tight_layout(pad=0)
def update(tp):  
    SkelMask = mainsks[tp,:,:]
    SkelMask = ndi.binary_dilation(SkelMask) #Dilation to make it more visible when ploting!
    SkelMask = np.ma.masked_where(SkelMask<1, SkelMask)
    plt.imshow(Maincrop[tp,:,:], interpolation='none', cmap = fluogreen)
    #plt.imshow(resmainmask[tp,:,:], interpolation='none', alpha = 0.2)
    plt.imshow(SkelMask, interpolation='none', cmap='Greys')
    
ani = FuncAnimation(fig, update, frames = time)  

writer = PillowWriter(fps=10)  
ani.save(GraphOutputFolder + r"\\" + filename[:-4]+  "mainaxism2.gif", writer=writer)  

print('we are saving minor axis of image')
fig = plt.figure(figsize=(10,2))
plt.axis('off')
fig.tight_layout(pad=0)
def update(tp):  
    SkelMask = minorsks[tp,:,:]
    SkelMask = ndi.binary_dilation(SkelMask) #Dilation to make it more visible when ploting!
    SkelMask = np.ma.masked_where(SkelMask<1, SkelMask)
    plt.imshow(Minorcrop[tp,:,:], interpolation='none', cmap = fluogreen)
    #plt.imshow(resminormask[tp,:,:], interpolation='none', alpha = 0.2)
    plt.imshow(SkelMask, interpolation='none', cmap='Greys')

ani = FuncAnimation(fig, update, frames = time)  

writer = PillowWriter(fps=10)  
ani.save(GraphOutputFolder + r"\\" + filename[:-4]+  r"minoraxism2.gif", writer=writer)  

#Saving graphs
print("Saving Graphs")
#Area graph:
plt.figure()
plt.title('Area variation')
plt.plot(realtime, ProjAreas, c="red")
plt.xlabel('Time (h)')
plt.ylabel('Area in pixels2')
plt.savefig(GraphOutputFolder + r"\\" + filename[:-4]+  'AreaVariationm2.png', bbox_inches='tight')
plt.show()
#Size graph

plt.figure()
plt.title('main and minor ax sizes')
plt.plot(realtime, MainAxSizes, c="red", label="Main axis")
plt.plot(realtime, MinorAxSizes, c="blue", label="Minor axis")
plt.xlabel('Time (h)')
plt.ylabel('Size in pixels')
plt.legend()
plt.savefig(GraphOutputFolder + r"\\" + filename[:-4]+  'AxSizeVariationm2.png', bbox_inches='tight')
plt.show()
#Curvature graph

plt.figure()
plt.title('main and minor ax mean curvatures')
plt.plot(realtime, mainmean, c="red", label="Main axis")
plt.plot(realtime, minormean, c="blue", label="Minor axis")
plt.xlabel('Time (h)')
plt.ylabel('Curvature')
plt.legend()
plt.savefig(GraphOutputFolder + r"\\" + filename[:-4]+ 'MeanCurvaturesm2.png', bbox_inches='tight')
plt.show()
#Exporting data
print("Exporting data")

allgood = input("Does it look good? Should we export the data?(y/n)")

if allgood == "y":
    #Creating output list for analysis!
    Outputlist = []
    Outputlist.append(SampleClass)
    Outputlist.append(AspectRatio)
    Outputlist.append(MainAxSizes)
    Outputlist.append(MinorAxSizes)
    Outputlist.append(MainSize)
    Outputlist.append(ProjAreas)
    Outputlist.append(mainmean)
    Outputlist.append(minormean)
    #Just in case adding original measured curvatures?
    Outputlist.append(maincurv)
    Outputlist.append(minorcurv)
    #In case we wanted to measure curvatures again with diff circle sizes:
    Outputlist.append(mainpoints)
    Outputlist.append(minorpoints)
    #Gathering a lot of data! some of it maybe unnecesary!
    
    #Great! Now we are ready to save it in specified folder and import to Pipeline 2
    
    #Save to output folder: 
        
    from datetime import date
    today = date.today()
    d1 = today.strftime("%Y_%m_%d") # This puts the date in format year, month day
    
    originalname = filename[:-4] #Taking filename without .tif to use it when saving data
    DataOutput = originalname + "_analyzed_" + d1 + ".csv"
    DataOut = join(DataOutputFolder, DataOutput)
    # We save the numpy array of surface area measurements for each timepoint 
    np.save(DataOut, Outputlist)
    print("Data has been saved to" + DataOut)
else:
    pass