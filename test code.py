#%%test code for github
#just adds and prints some numbers
import numpy as np

a = np.ones(10)
for i in np.arange(0,10):
    a[i] = a[i]+i**2

print(a)

