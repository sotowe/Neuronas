#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:38:50 2020

@author: jose
"""

# Load libraries
import numpy as np
import matplotlib.pyplot as plt


##------------------------------ 
#%%
r = np.linspace(0.05,1.4,100000)
eta1 = -1/(np.pi**2*r**2*4)-15*r+np.pi**2*r**2
eta2 = -1/(np.pi**2*r**2*4)-15*r+np.pi**2*r**2 - 3
v = 1/(np.pi*2*r)
plt.plot(eta1,r)
plt.plot(eta2,r)
plt.show()
plt.plot(eta1,v)
plt.plot(eta2,v)
plt.show()

v,r = np.meshgrid(np.linspace(-2.2,1.5,20),np.linspace(0.05,2,20))
etamedia = -5
J = 15
I = 0
rprima = 1/np.pi + 2*r*v
vprima = v**2+etamedia+J*r-np.pi**2*r**2+I
plt.quiver(r,v,rprima,vprima)
v = np.linspace(-2.2,1.5,200)
r = -1/(np.pi*v)
plt.plot(r,v)
r = np.linspace(-0.05,2,200)
v = -np.sqrt(-etamedia - J*r + np.pi**2*r**2)
plt.plot(r,v)
plt.xlim([0.05,2])
plt.ylim([-2.2,1.5])
plt.show()