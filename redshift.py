#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
H0 = 0.7*100 #km/s/Mpc
c = 300000000 #m/s
Dh = c/H0 
OmegaM = 0.3111
OmegaLambda = 0.6889
Omegak = 0

z = np.linspace(0,2,10000)
Dcint = 1/np.sqrt(OmegaM*(1+z)**3 + OmegaLambda)
Dc = Dh*scipy.integrate.cumtrapz(z,Dcint)
Da = Dc/(1+z[0:-1])
plt.plot(z[0:-1],Da/min(Da))
plt.xlabel("Redshift z")
plt.ylabel("Relative angular diameter distance")
# %%
