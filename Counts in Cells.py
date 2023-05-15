#%%
import matplotlib.pyplot as plt
import numpy as np
import math

#PROJECT CODE AVAILABLE AT https://github.com/quasur/GaBSe

#Common regions used
#format minRA maxRA minDec maxDec
bulletCellSize = 0.075
targetBulletRegion = 0.7
bulletmodulothingy =targetBulletRegion%bulletCellSize
bulletregionSize = targetBulletRegion-bulletmodulothingy#calculates scaling depending on target cellsize and region
bulletBoundary = np.array([100,110,-60,-5])
uniformBulletBoundary = np.array([103.57,103.57+bulletregionSize,-56.84,-56.84+bulletregionSize])

spteCellSize = 0.075
targetSpteRegion = 1.5
sptemodulothingy =targetSpteRegion%spteCellSize
spteregionSize = targetSpteRegion-sptemodulothingy
spteBoundary = np.array([71,76,-60,-55])
uniformSpteBoundary = np.array([73.5,73.5+spteregionSize,-58,-58+spteregionSize])

targetEGRegion = 0.5
elgordoCellSize = 0.075
EGmodulothingy =targetEGRegion%elgordoCellSize
EGregionSize = targetEGRegion-EGmodulothingy
elgordoBoundary = np.array([13.8,17.8,-50.5,-48.2])
uniformElgordoBoundary = np.array([14.732,14.732+EGregionSize,-49.194,-49.194+EGregionSize])


boundary = spteBoundary
uniformBoundary = uniformSpteBoundary
cellsize = spteCellSize
#Import data exported to .npy in datasetToNumpyConverter.py
with open("dataset.npy","rb") as f:
    datasetImport = np.load(f)
    
with open("reddata.npy","rb") as g: #import redmapper data
    redData = np.load(g)

#To center RA=0 in the plots
for i in range(np.size(datasetImport[0,:])):
    if datasetImport[0,i] > 160:
        datasetImport[0,i] = datasetImport[0,i] -360


#Function to cut out a region of the dataset based on min and max ra/dec coordinates. ra and dec coords supplied in arrays listed above
def regionMask(array,pos):
    ramin=pos[0]
    ramax=pos[1]
    decmin=pos[2]
    decmax=pos[3]
    array= array.T
    mask = np.logical_or(array[:,0]<ramin,array[:,0]>ramax)
    mask |= np.logical_or(array[:,1]<decmin,array[:,1]>decmax)
    retarray = array[~mask]
    return retarray.T

#Funtion to trim data to exclude points dimmer than the magnitude turnover point.
#These values in the maskX lines are calculated in colourGraphing.py
def magTurnoverCutoff(array):
    g = array[2]
    r = array[3]
    i = array[4]
    z = array[5]
    array=array.T

    maskG = np.logical_or(g>24.588131947517397,False)
    maskR = np.logical_or(r>24.103052185058594,False)
    maskI = np.logical_or(i>23.711337409973144,False)
    maskZ = np.logical_or(z>23.34798056793213,False)

    maskAll1 = np.logical_or(maskG,maskR)
    maskAll2 = np.logical_or(maskI,maskZ)
    maskAll = np.logical_or(maskAll1,maskAll2)

    gc = array[~maskG].T
    rc = array[~maskR].T
    ic = array[~maskI].T
    zc = array[~maskZ].T
    magsc = array[~maskAll].T

    return gc,rc,ic,zc,magsc

#Values higher than the turnouver point make up incomplete data as dim objects are harder to detect.
Gc,Rc,Ic,Zc,mag = magTurnoverCutoff(datasetImport)

#The largest of these cutoffs is in the I band. though they have similar quantities of data.
#print(np.size(Gc),np.size(Rc),np.size(Ic),np.size(Zc))
dataset = Ic

"""
plot entire data footprint"""

footprint = dataset
plt.rcParams.update({'font.size': 22})
plt.tight_layout()
plt.figure(figsize=(19,7))
plt.xlim(160,-30)
plt.ylim(-65,5)
plt.title("DES Data footprint")
plt.xlabel("Right Ascension (RA)/deg")
plt.ylabel("Declination (DEC)/deg")
plt.plot(footprint[0,:],footprint[1,:],'b.')
plt.savefig("graphs/LargeFootprint.png",format ="png")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
"""
plot smaller region"""

footprint = regionMask(dataset,boundary)
footprintcutout = regionMask(dataset, uniformBoundary)

minxgraph = min(footprint[0,:])
maxxgraph = max(footprint[0,:])
minygraph = min(footprint[1,:])
maxygraph = max(footprint[1,:])


plt.figure(figsize=((maxxgraph-minxgraph)*4,(maxygraph-minygraph)*4))
plt.rcParams.update({'font.size': 26})
plt.tight_layout()
plt.title("El Gordo Region")
plt.xlim(min(footprint[0,:]),max(footprint[0,:]))
plt.ylim(min(footprint[1,:]),max(footprint[1,:]))
plt.xlabel("Right Ascension/deg")
plt.ylabel("Declination/deg")
plt.plot(footprint[0,:],footprint[1,:],'b,')
plt.plot(footprintcutout[0,:],footprintcutout[1,:],'r,')
plt.savefig("graphs/El Gordo.png",format ="png")
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

#set region to run counts in cells on these set regions:
data = regionMask(dataset,boundary)
backgrounddata = regionMask(dataset,uniformBoundary)

def CIC(datacic,cellwidth,cellheight,fieldvar,fieldmean):
    #data is a NxM array where N contains x and y variables in columns 0 and 1 M is the number of points.
    #calculate region size and other useful values
    maxx = max(datacic[0,:])
    minx = min(datacic[0,:])
    maxy = max(datacic[1,:])
    miny = min(datacic[1,:])
    width = maxx-minx
    height = maxy-miny

    #determine the amount of cells needed - may have edge issues when bin size and width dont perfectly divide.
    binNumX = math.ceil(width/cellwidth)
    binnumY = math.ceil(height/cellheight)

    #translate data to origin to allow counts method to work
    datacic[0,:] = datacic[0,:] - minx
    datacic[1,:] = datacic[1,:] - miny

    #initialise array
    cellBins = np.zeros([binNumX,binnumY])

    #counts points in each bin
    for i in range(np.size(datacic[0,:])):
        xbin =math.floor(datacic[0,i]/cellwidth)
        ybin =math.floor(datacic[1,i]/cellheight)
        #print(xbin,ybin)
        cellBins[xbin,ybin] += 1
    
    #if there is no variance and mean provided, create one:
    if fieldvar == 0:
        fieldmean = np.sum(cellBins)/np.size(cellBins)        
        fieldvar = np.std(cellBins)

    significance = (cellBins - fieldmean)/fieldvar

    return significance,fieldvar,fieldmean

#sample refers to the smaller uniform sample within the larger field used to provide a good estimate of field mean and variance
sampleSig, sampleVar, sampleMean = CIC(backgrounddata,cellsize,cellsize,0,0)

fieldSig,fieldVar,fieldMean = CIC(data,cellsize,cellsize,sampleVar,sampleMean)

#if significance is above a certain value the cells that meet this have value 1 rather than 0
threshold = 3
clusters=fieldSig*0
clusters[fieldSig>threshold]=1 #sometimes my genius... its almost frightening
"""
#Plot image of significance of sample and field and identified clusters"""
fig,ax = plt.subplots(1,3,figsize=(30,10))
fig.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
plt.tight_layout()
smol =ax[0].imshow(sampleSig.T,origin="lower",cmap="afmhot")
ax[0].set_title("Significance in sample")
ax[0].set_ylabel("Cells/" + str(cellsize) + "deg")
ax[0].set_xlabel("Cells/" + str(cellsize) + "deg")
fig.colorbar(smol,ax=ax[0],shrink=0.8,pad=0.01)
ax[1].set_title("Significance of field")
big1 = ax[1].imshow(fieldSig.T,origin="lower",cmap="afmhot")
ax[1].set_ylabel("Cells/" + str(cellsize) + "deg")
ax[1].set_xlabel("Cells/" + str(cellsize) + "deg")
fig.colorbar(big1,ax=ax[1],shrink=0.8,pad=0.01)
ax[2].set_title("Candidate clusters")
big2 = ax[2].imshow(clusters.T,origin="lower",cmap="Greys")
ax[2].set_ylabel("Cells/" + str(cellsize) + "deg")
ax[2].set_xlabel("Cells/" + str(cellsize) + "deg")
plt.savefig("graphs/CICEG.png",format ="png")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

def clusterCutout(clusters,original,cellwidth,cellheight):
    
    maxx = max(original[0,:])
    minx = min(original[0,:])
    maxy = max(original[1,:])
    miny = min(original[1,:])

    clusterData = original[:,0:2]*0
    for i in range(np.size(clusters[:,0])):
        for j in range(np.size(clusters[0,:])):
            if clusters[i,j]==1:#for each cluster
                cellminx = i*cellwidth+minx
                cellminy = j*cellheight+miny
                cellmaxx = cellminx+cellwidth
                cellmaxy = cellminy+cellheight #calcuate the bounding region
                cellLoc = [cellminx,cellmaxx,cellminy,cellmaxy]
                currentClusterCell = regionMask(original,cellLoc)
                clusterData = np.append(clusterData,currentClusterCell,axis=1)
    clusterData = clusterData[:,2:]
    return clusterData

clusterData = clusterCutout(clusters,data,cellsize,cellsize)

#create a set of 3 random cells
randFakeClusters = np.zeros(np.shape(clusters))
randFakeClusters[20,18]=1 #random cell to compare to el gordo
randFakeClusterData = clusterCutout(randFakeClusters,data,cellsize,cellsize)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
"""
plot identified clusters"""
minxgraph = min(footprint[0,:])
maxxgraph = max(footprint[0,:])
minygraph = min(footprint[1,:])
maxygraph = max(footprint[1,:])
plt.figure(figsize=((maxxgraph-minxgraph)*4,(maxygraph-minygraph)*4))
plt.title("Identified clusers (red)")
plt.xlabel("Right Ascension/deg")
plt.ylabel("Declination/deg")
plt.xlim(0,max(data[0,:]))
plt.ylim(0,max(data[1,:]))
plt.plot(data[0,:],data[1,:],'b,')
plt.plot(clusterData[0,:],clusterData[1,:],'r,')
plt.savefig("graphs/EG Candidates.png",format ="png")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

clusterIZ = clusterData[4,:]-clusterData[5,:] #calculate colour
clusterIZ = clusterIZ[~(clusterIZ<-50)]#ensure the standard deviation isnt messed up by mags being at 99
randIZ = randFakeClusterData[4,:]-randFakeClusterData[5,:]
randIZ = randIZ[~(randIZ<-50)]
fieldIZ = data[4,:]-data[5,:]
mask = fieldIZ<-50
fieldIZ = fieldIZ[~mask]
print(np.std(clusterIZ))
print(np.std(fieldIZ))
print(np.std(randIZ))
"""
PLOT COLOUR MAGNITUDE DIAGRAM"""
plt.rcParams.update({'font.size': 26})
plt.figure(figsize=(16,16))
plt.plot(data[4,:],data[4,:]-data[5,:],'b,')#plot field colour data
plt.plot(clusterData[4,:],clusterData[4,:]-clusterData[5,:],'r.')#plot cluster colour data
plt.xlabel("I/mag")
plt.ylabel("I-Z/mag")
plt.title("Colour Magnitude diagram of candidate clusters")
plt.ylim(-3,2)
plt.xlim(15,24)
plt.savefig("graphs/SPTE Colour Mag.png",format="png")

redDataTrim = regionMask(redData,boundary) #trim redmapper data to the current regtion

clusterDEC = np.zeros(int(np.sum(clusters)))#initialise arrays
clusterRA = np.zeros(int(np.sum(clusters)))
clusterSig = np.zeros(int(np.sum(clusters)))
clusterMatch = np.zeros(int(np.sum(clusters)))

k=0
for i in range(np.size(clusters[:,0])):
    for j in range(np.size(clusters[0,:])):#for each cell
        if clusters[i,j]==1: #if there is a cluster present
            clusterRA[k]=boundary[0]+cellsize*(i+0.5)
            clusterDEC[k]=boundary[2]+cellsize*(j+0.5) #calculate center of cell
            clusterSig[k]=fieldSig[i,j]
            diffRA = np.abs(redDataTrim[0,:]-clusterRA[k])
            diffDEC = np.abs(redDataTrim[1,:]-clusterDEC[k])#calculates distance to current point to all redmapper points
            for l in range(np.size(diffRA)):
                if(diffRA[l]<=cellsize*1.5): #if this distance lies within the 9 bounding cells then it is counted as a match
                    if(diffDEC[l]<=cellsize*1.5):
                        clusterMatch[k]=1 #make match column =1 if match is success
            k=k+1
            
candidateData = np.vstack((clusterRA,clusterDEC,clusterSig,clusterMatch))




#plot redmapper data to candidate clusters. green is a match within 9 cells around a candidate, red is non match
plt.figure(figsize=(10,10))
plt.rcParams.update({'font.size': 26})
plt.xlim(71,76)
plt.ylim(-60,-55)
plt.xlabel("Right Ascension/deg")
plt.ylabel("Declination/deg")
plt.plot(redDataTrim[0,:],redDataTrim[1,:],'b.')
for i in range(np.size(candidateData[3,:])):
    if int(candidateData[3,i])==1:    
        plt.plot(candidateData[0,i],candidateData[1,i],'g.')
    elif int(candidateData[3,i])==0:
        plt.plot(candidateData[0,i],candidateData[1,i],'r.')
plt.savefig("graphs/redmapper comparison.png",format="png")
print(np.sum(candidateData[3,:])/np.size(candidateData[3,:]))


filename = "cells"+str(cellsize)+" candidates.npy"
with open(filename,"wb") as f:
    np.save(f,candidateData)
