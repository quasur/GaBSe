#%%adasd
import numpy as np
from astropy.io import fits

filename = 'sva1_gold_r1.0_catalog.fits'
# Open the given fits file
hdulist = fits.open(filename)
scidata = hdulist[1].data #turn fits into a managable table format
num = scidata.field(0) #cut array elements from the fits file
RA = scidata.field(1)
DEC = scidata.field(2)
cringe = scidata.field(3) #these 2 are parameters that say if the data could be bad or not
cringe2 = scidata.field(8)

data = np.vstack([RA,DEC,cringe,cringe2]) #create a numpy array of the desired data

def cringeMask(array): #Mask to reject points wich are classified as stars or dodgy
    array = array.T
    mask = np.logical_or(data[2]>1,data[3]>1) #if the point is a star or poor quality its rejected
    retarray = array[~mask]
    return retarray.T

expdata = cringeMask(data)

with open("desdata.npy","wb") as f:
    np.save(f,expdata)


#eee
