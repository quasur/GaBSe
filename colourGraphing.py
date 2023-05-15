#%%
#plots colour magnitude diagrams from DES SVA1 GOLD dataset
#https://github.com/quasur/GaBSe
import numpy as np
import matplotlib.pyplot as plt

with open("mag.npy","rb") as f:
    mags = np.load(f)

#if points have mag=99 to indicate they were undectected remove them from the array
def brightMask(array):
    array= array.T
    mask = np.logical_or(array[:,0]>95,array[:,1]>95)
    mask |= np.logical_or(array[:,2]>95,array[:,3]>95)
    retarray = array[~mask]
    return retarray.T


mags = brightMask(mags)
G = mags[0,:]
R = mags[1,:]
I = mags[2,:]
Z = mags[3,:]

#plot histograms of magnitude for each colour band
fig,ax = plt.subplots(2,2,figsize=(12,8))
fig.suptitle("\"Non-Galaxy\" data")

ax[0,0].set_title("G number histogram")
ax[0,0].set_yscale('log')
#ax[0,0].set_xlabel("Magnitude")
ax[0,0].set_ylabel("Number")
gd = ax[0,0].hist(G,1000,color="g")

ax[0,1].set_title("R number histogram")
#ax[0,1].set_xlabel("Magnitude")
ax[0,1].set_ylabel("Number")
ax[0,1].set_yscale('log')
rd = ax[0,1].hist(R,1000,color="orange")

ax[1,0].set_title("I number histogram")
ax[1,0].set_xlabel("Magnitude")
ax[1,0].set_ylabel("Number")
ax[1,0].set_yscale('log')
iid = ax[1,0].hist(I,1000,color="r")

ax[1,1].set_xlabel("Magnitude")
ax[1,1].set_ylabel("Number")
ax[1,1].set_title("Z number histogram")
ax[1,1].set_yscale('log')
zd =ax[1,1].hist(Z,1000,color="purple")


#%%
gda = gd[0]
gdb = gd[1]
index = int(np.where(gda==max(gda))[0])
print("max G",gdb[index])
rda = rd[0]
rdb = rd[1]
index = int(np.where(rda==max(rda))[0])
print("max R",rdb[index])
iida = iid[0]
iidb = iid[1]
index = int(np.where(iida==max(iida))[0])
print("max I",iidb[index])
zda = zd[0]
zdb = zd[1]
index = int(np.where(zda==max(zda))[0])
print("max Z",zdb[index])