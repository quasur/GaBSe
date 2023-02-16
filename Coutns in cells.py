#%%
import numpy as np
import matplotlib.pyplot as plt
import math

xrow = 8
yrow = 8

width = 1
dx = width/xrow
dy = width/yrow

cellcount = np.zeros([xrow,yrow])

with open("simpleset.npy","rb") as f:
    dataset = np.load(f)
testset = dataset[0,:]

for i in range(np.size(dataset[0,:])):
    xbin = math.floor(dataset[0,i]/dx)
    ybin = math.floor(dataset[1,i]/dy)
    cellcount[xbin,ybin] += 1

fieldVar = np.std(cellcount)**2
fieldMean = np.mean(cellcount)

cellVar = np.zeros([xrow,yrow])
clusterBool = np.zeros([xrow,yrow])
for i in range(xrow):
    for j in range(yrow):
        cellVar[i,j] = (cellcount[i,j]-fieldMean)/fieldVar
        if cellVar[i,j]>0:
            clusterBool[i,j] = 1



fig,ax = plt.subplots()
im = plt.imshow(cellcount,origin="lower")
for i in range(xrow):
    for j in range(yrow):
        text = ax.text(i,j,int(cellcount[i,j]),
                       ha="center", va="center", color="w")
plt.show()
plt.cla()
plt.imshow(clusterBool,origin="lower")


