#%%

#Opens .dat data from redMaPPer and converts it to .npy with relevant positional data
#https://github.com/quasur/GaBSe
import numpy as np
import pandas as pd

filename = 'sva1exp.dat'
# Open the given file
df = pd.read_csv(filename, sep="|", header=None)


RA = pd.Series.to_numpy(df[2])
DEC = pd.Series.to_numpy(df[3])

data = np.vstack((RA,DEC))

import matplotlib.pyplot as plt#plot footprint of reddata
plt.plot(data[0,:],data[1,:],'b,')

with open("reddata.npy","wb") as f:
    np.save(f,data)