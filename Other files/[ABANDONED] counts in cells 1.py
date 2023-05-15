#%% open des data
import numpy as np
import matplotlib.pyplot as plt
import math

#open data
with open("dataset.npy","rb") as f:
    dataset = np.load(f)

for i in range(np.size(dataset[0,:])):
    if dataset[0,i] > 160:
        dataset[0,i] = dataset[0,i] -360

        

def regionMask(array,ramin,ramax,decmin,decmax):
    array= array.T
    mask = np.logical_or(array[:,0]<ramin,array[:,0]>ramax)
    mask |= np.logical_or(array[:,1]<decmin,array[:,1]>decmax)
    retarray = array[~mask]
    return retarray.T

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

    gc = array[~maskG].T
    rc = array[~maskR].T
    ic = array[~maskI].T
    zc = array[~maskZ].T

    return gc,rc,ic,zc

Gc,Rc,Ic,Zc = magTurnoverCutoff(dataset)
#Ic is the largest
dataset = Ic

#%% plot footprint --------------------------------
plt.figure(figsize=(19,7))
plt.xlim(160,-30)
plt.ylim(-65,5)
plt.title("Data footprint")
plt.xlabel("RA/deg")
plt.ylabel("DEC/deg")
plt.plot(dataset[0,:],dataset[1,:],'b.')

#%%trim data to bullet----------------------------

bullet = regionMask(dataset,100,110,-60,-5)


minxg = min(bullet[0,:])
maxxg = max(bullet[0,:])
minyg = min(bullet[1,:])
maxyg = max(bullet[1,:])

#plot bullet cluster
plt.figure(figsize=((maxxg-minxg)*4,(maxyg-minyg)*4))
plt.title("Bullet cluster")
plt.plot(bullet[0,:],bullet[1,:],'b,')
plt.savefig("Bullet cluster.svg",format ="svg",dpi=1200)
plt.show()

# Bullet sample cutout-------------------------------------

#103.57
#104.28

#56.84
#56.13

bulletSample = regionMask(dataset,103.57,104.28,-56.84,-56.13)
plt.plot(bulletSample[0,:],bulletSample[1,:],'b,')
fieldVar = np.std(bulletSample)


#%% large block from SPT-E

bulk = regionMask(dataset,70,75,-60,-55)


minxg = min(bulk[0,:])
maxxg = max(bulk[0,:])
minyg = min(bulk[1,:])
maxyg = max(bulk[1,:])

#plot it
plt.figure(figsize=((maxxg-minxg)*4,(maxyg-minyg)*4))
plt.title("Region in SPT-E")
plt.plot(bulk[0,:],bulk[1,:],'b,')
plt.show()

# bulk C sample variance

sampleC = regionMask(dataset,73.67,73.67+1.25,-57.90,-57.90+1.33)
plt.plot(sampleC[0,:],sampleC[1,:],'b,')
fieldVar = np.std(sampleC)

#%%   COUNTS IN CELLS -----------

dataset = bullet

minx = min(dataset[0,:])
maxx = max(dataset[0,:])
miny = min(dataset[1,:])
maxy = max(dataset[1,:])

#allign with 0,0
dataset[0,:] = dataset[0,:]-minx
dataset[1,:] = dataset[1,:]-miny

width  = maxx-minx
height = maxy-miny


#loop to interate over various cell sizes
rowcounts = np.array([50])
for cellNum in rowcounts:

    xrows = cellNum
    yrows = int(cellNum/2)
    
    dx = width/xrows
    dy = height/yrows

    cellcount = np.zeros([xrows,yrows])
    
    sampleCArea = 1.25*1.33
    bulletSampleArea = (104.28-103.57)*(56.84-56.13)
    cellArea = dx*dy

    for i in range(np.size(dataset[0,:])):
        xbin = math.floor(dataset[0,i]/dx)-1
        ybin = math.floor(dataset[1,i]/dy)-1
        cellcount[xbin,ybin] += 1

    #fieldVar = np.std(cellcount)
    fieldMean = np.sum(cellcount)/np.size(cellcount)
    print(np.size(cellcount))
    #fieldMean = np.size(bulletSample[0,:])*cellArea/bulletSampleArea

    cellVar = np.zeros([xrows,yrows])
    clusterBool = np.zeros([xrows,yrows])
    for i in range(xrows):
        for j in range(yrows):
            cellVar[i,j] = (cellcount[i,j]-fieldMean)/fieldVar
            if cellVar[i,j]>3:
                clusterBool[i,j] = 1
    
    #plots an image with labels to show the cell count
    plt.rcParams.update({'font.size': 22})
    fig,ax = plt.subplots(2,1,figsize=(16,16))
    cic =ax[0].imshow(cellcount.T,origin="lower",cmap="afmhot")
    ax[0].set_title("Counts of cells")
    fig.colorbar(cic,ax=ax[:])
    
    plt.cla()
    ax[1].set_title("Clusters with cell number ")
    ax[1].imshow(clusterBool.T,origin="lower",cmap="Greys")
    plt.show()
    print("cell width: ",dx)
    print("cell height: ",dy)

#%% cut out candidate clusters----------------------------------------------------
clustercount = int(np.sum(clusterBool))
clusterpos =np.zeros([2,clustercount])
n=0
for i in range(np.size(clusterBool[:,0])):
    for j in range(np.size(clusterBool[0,:])):
        if clusterBool[i,j] != 0:
            clusterpos[0,n] = i
            clusterpos[1,n] = j
            n=n+1
clusterpos =clusterpos.T

clusterset = np.array([[],[],[],[],[],[]])
for i in range(np.size(clusterpos[:,0])):
    clusterminx = dx*clusterpos[i,0]+minx
    clustermaxx = dx*(clusterpos[i,0]+1)+minx

    clusterminy = dy*clusterpos[i,1]+miny
    clustermaxy = dy*(clusterpos[i,1]+1)+miny
    clusterRegion =regionMask(Ic,clusterminx,clustermaxx,clusterminy,clustermaxy) #set of points within a single cell of a potential cluster
    clusterset=np.append(clusterset,clusterRegion,axis=1)


bulletRegion =regionMask(Ic,100,110,-60,-5)
bulk = regionMask(Ic,70,75,-60,-55)
plt.figure(figsize=(16,16))
#plt.plot(bulk[0,:],bulk[1,:],'b,')
plt.plot(bullet[0,:],bullet[1,:],'b,')
plt.plot(clusterset[0,:],clusterset[1,:],'r,')

#%%

def brightMask(array):
    array= array.T
    mask = np.logical_or(array[:,0]>95,array[:,1]>95)
    mask |= np.logical_or(array[:,2]>95,array[:,3]>95)
    retarray = array[~mask]
    return retarray.T

clusterMags = brightMask(clusterRegion[1:-1,:])
regionMags = brightMask(bulletRegion[1:-1,:])

clusterGR = clusterMags[0,:]-clusterMags[1,:]
regionGR = regionMags[0,:]-regionMags[1,:]
plt.figure(figsize=(12,12))
plt.plot(regionGR,regionMags[1,:],'b,')
plt.plot(clusterGR,clusterMags[1,:],'r.')
plt.ylabel("G-R/relative magnitude")
plt.xlabel("R/relative magnitude")
plt.title("G-R graph for bullet cluster (red)")

#%%

magsCluster = clusterRegion[1:-1,:]
bulletCluster = bulletRegion[1:-1,:]

mags = brightMask(bulletCluster)

G = mags[0,:]
R = mags[1,:]
I = mags[2,:]
Z = mags[3,:]

fig,ax = plt.subplots(2,2,figsize=(16,12))
fig.suptitle("Bullet data")

ax[0,0].set_title("G number histogram")
ax[0,0].set_yscale('log')
#ax[0,0].set_xlabel("Magnitude")
ax[0,0].set_ylabel("Number")
gd = ax[0,0].hist(G,20,color="g")

ax[0,1].set_title("R number histogram")
#ax[0,1].set_xlabel("Magnitude")
ax[0,1].set_ylabel("Number")
ax[0,1].set_yscale('log')
rd = ax[0,1].hist(R,20,color="orange")

ax[1,0].set_title("I number histogram")
ax[1,0].set_xlabel("Magnitude")
ax[1,0].set_ylabel("Number")
ax[1,0].set_yscale('log')
iid = ax[1,0].hist(I,20,color="r")

ax[1,1].set_xlabel("Magnitude")
ax[1,1].set_ylabel("Number")
ax[1,1].set_title("Z number histogram")
ax[1,1].set_yscale('log')
zd =ax[1,1].hist(Z,20,color="purple")