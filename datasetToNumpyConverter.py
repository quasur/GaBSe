#%%adasd
import numpy as np
from astropy.io import fits

#converts .fits DES SVA1 GOLD data to .npy format for other files to use
#https://github.com/quasur/GaBSe

filename = 'sva1_gold_r1.0_catalog.fits'
# Open the given fits file
hdulist = fits.open(filename)
scidata = hdulist[1].data #turn fits into a managable table format
RA = scidata.field(1)
DEC = scidata.field(2)
type = scidata.field(3) #these 5 are parameters that say if the data could be bad or not
flagG = scidata.field(4)
flagR = scidata.field(5)
flagI = scidata.field(6)
flagZ = scidata.field(7)
G = scidata.field(9)
R = scidata.field(10)
I = scidata.field(11)
Z = scidata.field(12)

data = np.vstack([RA,DEC,type,flagG,flagR,flagI,flagZ]) #create a numpy array of the desired data
magdata = np.vstack([G,R,I,Z,type,flagG,flagG,flagI,flagZ])

def cringeMask(array): #Mask to reject points wich are classified as non galaxies or dodgy
    array = array.T
    mask = np.logical_or(data[-5]!=1,data[-4]>0)
    mask |= np.logical_or(data[-3]>0,data[-2]>0)
    mask |= np.logical_or(data[-1]>0,False) #if the point is a star or poor quality its rejected
    retarray = array[~mask]
    return retarray.T

expdata = cringeMask(data)
magdata = cringeMask(magdata)#remove poor data from set
expdata = expdata[0:2,:]
magdata = magdata[0:4,:]

with open("desdata.npy","wb") as f: #export pos data
    np.save(f,expdata)

with open("mag.npy","wb") as g: #export colour data
    np.save(g,magdata)

dataset = np.vstack([expdata,magdata])

with open("dataset.npy","wb") as h: #combination of colour and positional data
    np.save(h,dataset)


#eee
