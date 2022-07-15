# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 15:15:45 2021

@author: ara
"""

##New pipeline for the analysis of deformation of human embryonic bodies
#Objectives, use mCherry channel to segment and measure changes in:
#Area, eccentriciy, perimeter and curvature of the perimeter.

  
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from skimage import filters, measure, morphology
import math
from tifffile import imwrite, TiffFile
#import tifffile
#import seaborn as sns

#First let's import image and define path
folder = input("Please introduce name of folder")
filename = input("Please introduce name of file you want to analyze, remember to add .tiff")
path = r'C:\Users\ara\Documents\Python Scripts\Embryonic body constriction\Samples'
from os.path import join

filepath = join(path, folder, filename)
I = io.imread(filepath)

print("file has been loaded")
tps = I.shape[0] #File shape is (timepoints, z, x, y, channels) 

red = I[:,:,:,:,1] #mCherry channel, we pick by default all slices to project
with TiffFile(filepath) as tif:
    image = tif.asarray()
    imagej_metadata = tif.imagej_metadata
    
xyres = imagej_metadata['spacing']

#Stimulation area parameters:in pixels, 0 is position of start of  square and 1 is end position.
#disctionary with all the stimulations for all the samples: 
stimareas = {"Sample 1" : [65, 272, 137, 373], 
             "Sample 2" : [165, 58, 244, 122], 
             "Sample 3" : [209, 378, 292, 456],
             "Sample 4" : [380, 152, 449, 214], 
             "Sample 5" : [208, 57, 268, 124], 
             "Sample 6": [386, 171, 460, 237],
             "Sample 7": [218, 390, 287, 455],
             "Sample 8": [200, 395, 312, 464]}
#They are ordered as y0, x0, y1, x1. For the first and second corner of the bounding box of the stim area

stimy0, stimx0, stimy1, stimx1 = stimareas[folder]
stimcenty =(stimy0 + stimy1)/2
stimcentx = (stimx0 + stimx1)/2
#By averaging each, we get the central position of the stimulation. 

#Step 1, get nice segmentation of red channel to then extract measurements

#From old adapted function

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
        blur = filters.gaussian(Maxproject[t,:,:], sigma = 5)
        otsu = filters.threshold_otsu(blur)
        mask = blur > otsu * 0.75
        mask = ndi.binary_opening(mask)
        mask = ndi.binary_closing(mask, iterations=3)
        mask = ndi.binary_fill_holes(mask)
        mask = ndi.binary_erosion(mask, iterations = 7)
        masks[t,:,:] = mask
    #Now that we have the mask, will have to select biggest one. To do so we label and get biggest area 
        labelled_all[t,:,:], nr = ndi.label(masks[t, :,:])
        props = measure.regionprops(labelled_all[t,:,:])
        area_t = np.array([prop['area'] for prop in props])
        if t == 0: #30thMarch2020 Added this to track and keep always the colony that is biggest at t=0
            biggest = max(enumerate(area_t), key=lambda x : x[1])[0] + 1
            indexs[t] = biggest
            biggestMask[t,:,:] = labelled_all[t,:,:] == biggest
            biggestMask[t,:,:] = ndi.binary_dilation(biggestMask[t,:,:] , iterations = 5)
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
            biggestMask[t,:,:] = ndi.binary_dilation(biggestMask[t,:,:] , iterations = 5)
    return Maxproject, biggestMask


print('Doing projection')
projections, mainmasks = ProjectAndSegment(red)
print('Projection done!')
#plt.imshow(mainmasks[1])
#plt.imshow(mainmasks[6]) #It looks great!

#It seems like it worked on the first trial, let's see what 
#measurements can we get with measure region props. 
area = []
ecc = []
per = []
time = [] #To plot real time
diameter = []
for t in range(tps):
    props = measure.regionprops(mainmasks[t])
    area.append(props[0].area*xyres**2)
    ecc.append(props[0].eccentricity)
    per.append(props[0].perimeter*xyres)
    time.append(t*5) #In minutes!
    diameter.append(props[0].major_axis_length*xyres)

#plt.figure()
#plt.plot(time, per)
#plt.show()

#Ok! There we have all the very basic measurements. Now what about curvature? 

#We will measure curvature from perimeter: 

perim = np.zeros((tps,512,512))
for t in range(tps):
    temp = ndi.binary_dilation(mainmasks[t])
    temp = temp - mainmasks[t]
    perim[t,:,:] = temp

#plt.imshow(perim[-4])  
#Alright we have perimeters now let's go counter clock

#My old macro function to put points in order: 

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
        
        nbnumber=np.zeros((0, 4)) #Will be [time, x, y, number of neighbours] index not needed because is kept on the first for iteration
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
        #Now we have the list with number of neighbours,
        #CHANGE HERE: Starting point will directly be the on most to the right. 
        majorx = np.amax(nbnumber[:,1])
        maxright = np.where(nbnumber[:,1] == majorx)
        if len(maxright[0]) == 1:
            firstpoint = nbnumber[np.where(nbnumber[:,1] == majorx),:]
            firstpoint = firstpoint[0][0]
        else:
            firstpoint = nbnumber[maxright[0][0],:]
        #print("first point is")
        #print(firstpoint)
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
                            neighbourscheck=np.concatenate((upleft[0], up[0], upright[0], left[0], right[0], downleft[0], down[0], downright[0]), axis=None) #it should always choose upwards first
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



#If that works properly, the list of points will go from 0 to 2pi. 
print('Ordering points, getting curvature')
mainpoints = PointsInOrder(perim)

#Checking that order is correct: 

#Forth row must be ascending 
for t in range(tps):
    where = np.where(mainpoints[:,0] == t)
    where = where[0]
    #print(where)
    yvalues = mainpoints[where,3]
    if yvalues[0] > yvalues[80]: #Question is Y ascending? if no
        mainpoints[where,:] = mainpoints[where[::-1],:] #Reverse those rows
        print("timepoint " + str(t) + " reverted")
    else:
        pass
#Old curvature functions

def RadToDeg(a):
    b=a*180/math.pi
    if b< 0:
        b=360 + b
    else:
        pass
    return b
   
def FindCurvature(points, Csize, Csep, xyres):
    Circlesize = Csize #Distance between the points of the circle
    Circlesep = Csep #How often we pick points for circles, with 1 we pick every pixel for a circle with its following 2 points. 
    frames = np.amax(points[:,0], axis=0) +1
    Curvatures = []
    for t in range(int(frames)):
        Pointstemp = points[np.where(points[:,0] == t)]
        CurvaturesTemp= [] #Will store Curvature
        Pointstemp2 = np.vstack((Pointstemp, Pointstemp)) # This way it we can iterate through the points at the beggining again
        for i in range(int((Pointstemp.shape[0])/Circlesep)):  #This way it will stop at the first point of the last circle
            cs=i*Circlesep  #Circle starts at
            pointA=np.array([Pointstemp2[cs-Circlesize,2], Pointstemp2[cs-Circlesize,3]])
            pointB=np.array([Pointstemp2[cs,2], Pointstemp2[cs,3]])
            pointC=np.array([Pointstemp2[cs+Circlesize,2], Pointstemp2[cs+Circlesize,3]])        
            sidea=math.sqrt(((pointA[0]-pointB[0])**2)+((pointA[1]-pointB[1])**2))
            sideb=math.sqrt(((pointB[0]-pointC[0])**2)+((pointB[1]-pointC[1])**2))
            sidec=math.sqrt(((pointC[0]-pointA[0])**2)+((pointC[1]-pointA[1])**2))
            s=(sidea+sideb+sidec)/2
            if s*(s-sidea)*(s-sideb)*(s-sidec) <= 0: #This happens when the three dots are in a line, which makes curvature == 0
                Curv = 0
            else:
                K=math.sqrt(s*(s-sidea)*(s-sideb)*(s-sidec))
                Curv=4*K/(sidea*sideb*sidec*xyres)#added resolution so curvature is in 1/um
            #Adding positive and negative 3rd April2020
            AB = [pointA[0]-pointB[0], pointA[1]-pointB[1]]
            BC = [pointB[0]-pointC[0], pointB[1]-pointC[1]]
            angleABC = math.atan2(BC[1], BC[0]) - math.atan2(AB[1], AB[0])
            angleABC = RadToDeg(angleABC)
            if angleABC >= 180: # if right turn, negative curvature
                Curv = -Curv
            else:
                pass
            CurvaturesTemp.append(Curv)
        Curvatures.append(CurvaturesTemp)
    return Curvatures
#How to define which circle size to pick?
#Circle size really impacts a lot, I should relate it to radius to see
#How much it aproaches the real curvature (Should oscilate around it)
#Since we know diameter of the sample
radius0 = diameter[0]/2
#in this case of 135
#How to have the appropiate resolution? From the radius measured, we expect a perimeter:
exper = 2*math.pi*radius0
#In this case of 848 microns.
# We want to measure curvatures of circles from, let's say 50 to 2000 (two orders of magnitude)
#Then the lowest curvature we want to be able to measure is the one of a circle of 50um radius
# which would have a perimeter of 62um, in pixels: 62/xyres
#Since circle size must be the half of it, because it represents the distance between 2 of the 3 points:
circlesize = int(400/2*xyres)
#This circle size, of approximately 100 pixels, will neglect small very small circles (high curvatures caused by artifacts) by measuring points that are further away of their smaller perimeter
#And will be big enough to accurately measure big circles(small curvatures) accurately, based on the tests I did previously. 
maincurv = FindCurvature(mainpoints, circlesize, 1, xyres)


#To plot curvature we will need a centroid for every timepoint. 
def getCenters(maskarray):
    #Mask array is a 3D array, of 2D masks in time from which we will get centroid for each
    #dims are [t,y,x]
    centroids = np.zeros([0,2])
    time, y, x = maskarray.shape
    for t in range(time):
        props = measure.regionprops(mainmasks[t])
        centroidy, centroidx = props[0].centroid
        centroids = np.vstack((centroids,[centroidy,centroidx]))
    return centroids

allcentroids = getCenters(mainmasks)
# Now that we have a list of centroids, we can measure the angle for each curvature measured on every point
#We will be considering the curvature measured on the central point of the osculating circle.

def getAngle(origin, point):
    #Origin and point are coordinates in y, x. 
    #Angle will be calculated based on origin provided.#
    angle = math.atan2(point[1]-origin[1],point[0]-origin[0]) # atan function is for an x, y point
    decimalangle = RadToDeg(angle)
    return decimalangle

curvsangles = []
for t in range(tps):
    currentcang = np.zeros((0 ,4))
    currentpoints = mainpoints[np.where(mainpoints[:,0] == t)]
    for i in range(currentpoints.shape[0]):
        tempangle = getAngle(allcentroids[t,:], currentpoints[i,2:4]) #for some reason, centroid coming as x, y??
        currentcang = np.vstack((currentcang, [tempangle, maincurv[t][i], currentpoints[i,2], currentpoints[i,3]]))
    curvsangles.append(currentcang)



#I would like to paint the curvature with the location of the point to validate the measurement
#will include, the image, the mask, the curvature, and the measured angle.
#To "paint" curvature not only on a single pixel, we will have to apply a mask on that location
def circleMask(radius, margin = 5):
    #Radius in pixels
    imsize = radius*2+margin*2
    mask = np.zeros([imsize,imsize])
    for x in range(imsize):
        for y in range(imsize):
            if math.sqrt((x-radius-margin)**2 + (y-radius-margin)**2)< radius:
                mask[y,x]=1
    return mask

def paint(array, kernel, location, value):
    #array is the array to paint on (np)
    #kernel is a binary numpy array with the shape to be painted
    #location is a 2 element list with y,x location for the center of the shape
    # value is the value to be painted, int or float
    ay, ax = array.shape
    ky, kx = kernel.shape
    krows, kcols = np.where(kernel == 1)
    paintrows = krows +location[0] - ky//2
    paintcols = kcols + location[1] - kx//2
    #Those are the coordinates to paint:
    for y in paintrows:
        for x in paintcols:
            if x>=ax or y >= ay:
                continue
            array[y,x] =  value
    #!!! It will modify the original array!
    return
        
circ = circleMask(4)
#From coords we can place a circle with that center on the image.
print("making images with curvature data")
angleimage = np.zeros_like(mainmasks)
curvatureimage = np.zeros_like(mainmasks, dtype = "float16")
for t in range(tps):
    print(t)
    for i in range(len(curvsangles[t])):
        paint(angleimage[t,:,:], circ, [int(curvsangles[t][i,2]), int(curvsangles[t][i,3])], curvsangles[t][i,0])
        paint(curvatureimage[t,:,:], circ, [int(curvsangles[t][i,2]), int(curvsangles[t][i,3])], curvsangles[t][i,1])

#Now we want to save the image:
tps, y, x = mainmasks.shape
finalimage = np.zeros((tps, 4, y, x))
finalimage[:,0,:,:] = projections
finalimage[:,1,:,:] = mainmasks
finalimage[:,2,:,:] = angleimage
finalimage[:,3,:,:] = curvatureimage
finalimage= finalimage.astype("float32")
imwrite(filepath[:-4]+ r"curvatureimage.tif", finalimage, imagej = True, metadata={'axes': 'TCYX'})

#This outputs a matrix of angle - curvature, 

#We manually input the data of where the stimulation area was at the beginning.
#From this we calculate the angle of the center of the stimulation, and consider this 0degrees
# and we will plot from -180 - 0 to 180. Let's get the angle:

stimangle = getAngle(allcentroids[0,:], [stimcenty, stimcentx])
#We just need to substract and make everything over 180 go to the other side:
for t in range(tps):
    for i in range(len(curvsangles[t][:,0])):
        curvsangles[t][i,0] = curvsangles[t][i,0] - stimangle
        if curvsangles[t][i,0] > 180:
            curvsangles[t][i,0] = curvsangles[t][i,0] - 360
        elif curvsangles[t][i,0] < -180:
            curvsangles[t][i,0] = curvsangles[t][i,0] + 360

t= 10 
figure, ax =plt.subplots()
ax.plot(curvsangles[t][:,0],curvsangles[t][:,1])
ax.set_ylabel("curvature (1/um)")
#maxlim = np.amax(curvsangles[t][:,1]) + 0.1*np.amax(curvsangles[t][:,1])
#minlim = np.amin(curvsangles[t][:,1]) - 0.1*np.amin(curvsangles[t][:,1])
#ax.set_ylim(minlim, maxlim)
#ax2 = ax.twinx()
#ax2.set_ylim(minlim, maxlim)
#ax2.set_ylabel("Osculating circle radius (um)")
plt.show()


# Guess we can manually group every 1 degree? 
curvsandangles = np.zeros((tps, 360, 2))
for t in range(tps):
    for i in range(360):
        low = i -180
        high = i+1 -180
        pickedcvs = []
        for a in range(curvsangles[t].shape[0]):
            if curvsangles[t][a,0] > low and curvsangles[t][a,0] < high:
                pickedcvs.append(curvsangles[t][a,1])
            
        avcurvature = np.mean(pickedcvs)
        curvsandangles[t,i,0] = (low+high)/2
        curvsandangles[t,i,1] = avcurvature

t1= 2
t2 = 14
t3 = -1
figure, ax =plt.subplots(figsize = (10,6))
ax.plot(curvsandangles[t1,:,0],curvsandangles[t1,:,1], color = "black", alpha = 0.8, label = "Before stimulation")
ax.plot(curvsandangles[t2,:,0],curvsandangles[t2,:,1], color = "blue", alpha = 0.8, label = "1h stimulation")
ax.plot(curvsandangles[t3,:,0],curvsandangles[t3,:,1], color = "red", alpha = 0.8, label = "15 min post stimulation")
ax.set_ylabel("curvature (1/um)")
ax.legend(loc = "lower left")
#maxlim = np.amax(curvsandangles[t2][:,1]) + 0.1*np.amax(curvsandangles[t2][:,1])
#minlim = np.amin(curvsandangles[t2][:,1]) - 0.1*np.amin(curvsandangles[t2][:,1])
#ax.set_ylim(minlim, maxlim)
#ax2 = ax.twinx()
#ax2.set_ylim(1/minlim,1/maxlim)
#ax2.set_ylabel("Osculating circle radius (um)")
plt.savefig(filepath[:-4] + 'maincurvatures.png', bbox_inches='tight')
plt.show()


#Gotta check that the angle is ok! to reorient the thing: 
#Display as circular plot,
#t = 3
#fig = plt.figure()
#ax = fig.add_subplot(projection='polar')
#c = ax.plot(curvsandangles[t][:,0]*math.pi/180,curvsandangles[t][:,1], color= 'blue')


#Anyway we can get all measurements and think of calculations in the next stage

#Let's add mean curvature so it's already calculated: 
meancvs = []
for t in range(tps):
    meancvs.append(np.mean(maincurv[t]))
    
#print(meancvs)
#Alright time to print report and store data!
plt.rcParams.update({'font.size': 16})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

plt.figure(figsize = (10, 16), dpi=200)
#plt.suptitle('Distribution of radial component in all vectors')
ax5 = plt.subplot(3,2,1)
plt.imshow(mainmasks[0])
ax6 = plt.subplot(3,2,2)
plt.imshow(perim[0])
ax1 = plt.subplot(3,2,3)
plt.title('Area')
plt.axvspan(15, 70, color = '#78e3d8', alpha = 0.3)
plt.plot(time, area, c='k')
plt.ylabel('area (pixels**2)')
plt.xlabel('time (min)')
#plt.xlabel('Distance to centre(um)')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2 = plt.subplot(3,2,4)
plt.title('Eccentricity')
plt.axvspan(15, 70, color = '#78e3d8', alpha = 0.3)
plt.plot(time, ecc, c='k')
plt.ylabel('eccentricity')
plt.xlabel('time (min)')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3 = plt.subplot(3,2,5)
plt.title('perimeter')
plt.axvspan(15, 70, color = '#78e3d8', alpha = 0.3)
plt.plot(time, per, c='k')
#plt.ylabel('velocity (um/min)')
plt.ylabel('perimeter (pixels)')
plt.xlabel('Time (min)')
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax4 = plt.subplot(3,2,6)
plt.title('average curvatures')
plt.axvspan(15, 70, color = '#78e3d8', alpha = 0.3)
plt.plot(time, meancvs, c='k')
plt.xlabel('Time (min)')
plt.ylabel('Curvature')
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
#plt.legend()
plt.savefig(filepath[:-4] + 'report.png', bbox_inches='tight')
plt.show()

from matplotlib.animation import FuncAnimation, PillowWriter  

#To save perimeter evolution
print('we are saving perimeter evolution')
fig = plt.figure(figsize=(10,10))
plt.axis('off')
fig.tight_layout(pad=0)
def update(tp):  
    plt.imshow(curvatureimage[tp,:,:], interpolation='none')
    
ani = FuncAnimation(fig, update, frames = tps)  

writer = PillowWriter(fps=10)  
ani.save(filepath[:-4] + "projectionsm2.gif", writer=writer) 



print("Exporting data")

allgood = input("Does it look good? Should we export the data?(y/n)")

if allgood == "y":
    #Creating output list for analysis!
    Outputlist = []
    Outputlist.append(area)
    Outputlist.append(ecc)
    Outputlist.append(per)
    #Just in case adding original measured curvatures?
    Outputlist.append(maincurv)
    Outputlist.append(meancvs)
    Outputlist.append(diameter)
    #In case we wanted to measure curvatures again with diff circle sizes:

    #Save to output folder: 
        
    from datetime import date
    today = date.today()
    d1 = today.strftime("%Y_%m_%d") # This puts the date in format year, month day
    
    originalname = filename[:-4] #Taking filename without .tif to use it when saving data
    DataOutputFolder = r"C:\Users\ara\Documents\Python Scripts\Embryonic body constriction\Output Data"
    DataOutput = originalname + "_analyzed_" + d1 + ".csv"
    DataOut = join(DataOutputFolder, DataOutput)
    # We save the numpy array of surface area measurements for each timepoint 
    np.save(DataOut, Outputlist)
    DataOutput2 = originalname + "curvatures_analyzed_" + d1 + ".csv"
    DataOut2 = join(DataOutputFolder, DataOutput2)
    np.save(DataOut2, curvsandangles)
    print("Data has been saved to" + DataOut)
else:
    pass
