#%% Generates samples to test counts in cells on
#https://github.com/quasur/GaBSe

import numpy as np
import matplotlib.pyplot as plt

#generating a random population with a region of high density
seed=1233
np.random.seed(seed)

sourceAu = np.vstack((np.random.normal(loc=0.25,scale=0.02,size=200) , np.random.normal(loc=0.25,scale=0.02,size=200)))
sourceBu = np.vstack((np.random.normal(loc=0.75,scale=0.08,size=200) , np.random.normal(loc=0.25,scale=0.08,size=200)))
sourceCu = np.vstack((np.random.normal(loc=0.25,scale=0.04,size=200) , np.random.normal(loc=0.75,scale=0.04,size=200)))
sourceDu = np.vstack((np.random.normal(loc=0.75,scale=0.10,size=200) , np.random.normal(loc=0.75,scale=0.10,size=200)))

def setMask(array): #Funcion I found to simply remove values outside of the unit square
    array= array.T
    xmin,ymin = 0,0
    xmax,ymax = 1,1
    mask = np.logical_or(array[:,0]<xmin,array[:,0]>xmax)
    mask |= np.logical_or(array[:,1]<ymin,array[:,1]>ymax)
    retarray = array[~mask]
    return retarray.T


sourceA = setMask(sourceAu)
sourceB = setMask(sourceBu)
sourceC = setMask(sourceCu)
sourceD = setMask(sourceDu)

lic = 1200
#if you want a more obvious cluster, increase hic or decrease the 500 in lic
lowIntensity = setMask(np.array([np.random.rand(lic),np.random.rand(lic)])) #generate a random set of points
#plot sources, change the "r." to "b." if you want to see the 2 types clearer
plt.figure(figsize=(6,6))
plt.plot(lowIntensity[0,:],lowIntensity[1,:],color="red",linestyle="none",marker=".")
plt.plot(sourceA[0,:],sourceA[1,:],'y.')
plt.plot(sourceB[0,:],sourceB[1,:],'g.')
plt.plot(sourceC[0,:],sourceC[1,:],'b.')
plt.plot(sourceD[0,:],sourceD[1,:],color="purple",linestyle="none",marker=".")
plt.title("Generated Testing Samples")
plt.savefig("graphs/Fake Samples.png",format="png")



dataset = np.hstack((sourceA,sourceB,sourceC,sourceD,lowIntensity))
with open("gaussianset.npy","wb") as f:
    np.save(f,dataset)


