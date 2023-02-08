#%% Counts in cells algorithm for finding high regions of surface density

import numpy as np
import matplotlib.pyplot as plt

#generating a random population with a region of high density
np.random.seed(1234)

hic = 1500 #high intensity count
highIntensity = np.array([np.random.rand(hic),np.random.rand(hic)])

#circle of r = 0.2 centred at 0.4,0.4 = =>
# (x-0.4)^2 + (y-0.4)^2 = 0.2^2
inCircle = []#initialise array

#checks if each point is within the circle and discards ones that arent
for i in range(hic):
    circleCheck = ((highIntensity[0,i]-0.4)**2+(highIntensity[1,i]-0.4)**2)
    if circleCheck < 0.04:
        inCircle.append(highIntensity[:,i])
clusterSource = np.transpose(np.array(inCircle))
#likely a faster way to do this but its not worth it at this low sample size


lic = 500-np.size(clusterSource[0,:])  #low intensity count, using 500 to keep a round number
#if you want a more obvious cluster, increase hic or decrease the 500 in lic
lowIntensity = np.array([np.random.rand(lic),np.random.rand(lic)]) #generate a random set of points

#plot sources, change the "r." to "b." if you want to see the 2 types clearer
plt.plot(clusterSource[0,:],clusterSource[1,:],'r.')
plt.plot(lowIntensity[0,:],lowIntensity[1,:],'r.')
plt.show()

