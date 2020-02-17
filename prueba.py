#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 09:54:46 2020

@author: jose
"""

# Load libraries
import numpy as np
import matplotlib.pyplot as plt


##------------------------------ 
#%%
#%%

# Estudiamos los puntos estables para un J determinado en función de etamedia

J = 15
# Primero r en función de eta. 
r = np.linspace(0.05,1.6,100000)

# Sin intensidad de corriente
eta1 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2

# Con intensidad de corriente
eta2 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2 - 3

# Ahora v en función de eta
v = -1/(np.pi*2*r)
a = np.loadtxt('archivostxt/meanaverageV_J15_ro0_Vo-3.txt')
b = np.loadtxt('archivostxt/meanaverageV_J15_ro1_Vo0.txt')
# Ploteamos las funciones
plt.plot(eta1,r)
plt.xlabel("eta")
plt.ylabel("r")
plt.title("J=15")
plt.xlim([-11,0])
plt.savefig('r_vs_eta_J15.png') 
plt.show()
plt.plot(eta1,v)
plt.plot(a[:,0],a[:,1], 'bo')
plt.plot(b[:,0],b[:,1], 'bo')
plt.xlabel("eta")
plt.ylabel("v")
plt.title("J=15")
plt.xlim([-11,0])
plt.savefig('v_vs_eta_J15.png') 
plt.show()
