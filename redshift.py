#%%
#Calculates the expected angular size of galaxy clusters dependent on redshift.
#https://github.com/quasur/GaBSe
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
H0 = 0.7*100 #km/s/Mpc
c = 300000 #km/s
Dh = c/H0 
OmegaM = 0.3111
OmegaLambda = 0.6889
Omegak = 0

z = np.linspace(0,3,10000)
Dcint = 1/np.sqrt(OmegaM*(1+z)**3 + OmegaLambda)
Dc = Dh*scipy.integrate.cumtrapz(Dcint,z)
Da = Dc/(1+z[0:-1])*(np.pi/180)
Deg = 1/Da

#plotting functions
plt.rcParams.update({'font.size': 22})
fig,ax = plt.subplots(1,3,figsize=(27,9))
ax[2].plot(z[0:-1],Deg)
ax[2].set_ylim(0,0.1)
ax[2].set_xlim(0,1)
ax[2].set_xlabel("Redshift z")
ax[2].set_ylabel("Angular size/degree")
fig.suptitle("Angular size of typical galaxy cluster")
ax[1].plot(z[0:-1],Deg)
ax[1].set_xlim(0,3)
ax[1].set_ylim(0,0.3)
ax[1].set_xlabel("Redshift z")
ax[1].set_ylabel("Angular size/degree")
ax[0].plot(z[0:-1],Da)
ax[0].set_xlim(0,3)
ax[0].set_xlabel("Redshift z")
ax[0].set_ylabel("Angular Diameter Distance/[Mpc/degree]")
fig.subplots_adjust(left=0.04, right=0.96, bottom=0.1, top=0.9)
fig.savefig("graphs/Redshift.png",format ="png")