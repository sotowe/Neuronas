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

# Estudiamos los puntos estables para un J determinado en función de etamedia

J =5
# Primero r en función de eta. 
r = np.linspace(0.04,0.75,100000)

# Sin intensidad de corriente
eta1 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2

# Con intensidad de corriente
eta2 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2 - 3

# Ahora v en función de eta
v = -1/(np.pi*2*r)

# cargamos ficheros de texto con los datos experimentales
a1 = np.loadtxt('archivostxt/meanrate_J05_ro1_Vo0.txt')
# a2 = np.loadtxt('archivostxt/meanrate_J05_ro1_Vo0_2.txt')
a3 = np.loadtxt('archivostxt/meanrate_J05_ro0_Vo-3.txt')
# a4 = np.loadtxt('archivostxt/meanrate_J055_ro0_Vo-3_2.txt')
# =============================================================================
b1 = np.loadtxt('archivostxt/meanaverageV_J05_ro1_Vo0.txt')
# b2 = np.loadtxt('archivostxt/meanaverageV_J05_ro1_Vo0_2.txt')
b3 = np.loadtxt('archivostxt/meanaverageV_J05_ro0_Vo-3.txt')
# b4 = np.loadtxt('archivostxt/meanaverageV_J05_ro0_Vo-3_2.txt')
# =============================================================================
# Ploteamos las funciones
# Para r
plt.plot(eta1,r)
plt.errorbar(a1[:,0],a1[:,1], a1[:,2], fmt = 'bo')
# plt.errorbar(a2[:,0],a2[:,1],a2[:,2],fmt = 'bo')
plt.errorbar(a3[:,0],a3[:,1],a3[:,2],fmt = 'ro')
# plt.errorbar(a4[:,0],a4[:,1],a4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("r")
plt.title("J=5")
plt.xlim([-11,0])
plt.savefig('r_vs_eta_J5.png') 
plt.show()

#Para v
plt.plot(eta1,v)
plt.errorbar(b1[:,0],b1[:,1], b1[:,2], fmt = 'bo')
# plt.errorbar(b2[:,0],b2[:,1], b2[:,2], fmt = 'bo')
plt.errorbar(b3[:,0],b3[:,1], b3[:,2], fmt = 'ro')
# plt.errorbar(b4[:,0],b4[:,1], b4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("v")
plt.title("J=5")
plt.xlim([-11,0])
plt.savefig('v_vs_eta_J5.png') 
plt.show()


#%%

# Estudiamos los puntos estables para un J determinado en función de etamedia

J =10
# Primero r en función de eta. 
r = np.linspace(0.04,1.2,100000)

# Sin intensidad de corriente
eta1 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2

# Con intensidad de corriente
eta2 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2 - 3

# Ahora v en función de eta
v = -1/(np.pi*2*r)

# cargamos ficheros de texto con los datos experimentales
a1 = np.loadtxt('archivostxt/meanrate_J10_ro1_Vo0.txt')
a2 = np.loadtxt('archivostxt/meanrate_J10_ro1_Vo0_2.txt')
a3 = np.loadtxt('archivostxt/meanrate_J10_ro0_Vo-3.txt')
a4 = np.loadtxt('archivostxt/meanrate_J10_ro0_Vo-3_2.txt')
# =============================================================================
b1 = np.loadtxt('archivostxt/meanaverageV_J10_ro1_Vo0.txt')
b2 = np.loadtxt('archivostxt/meanaverageV_J10_ro1_Vo0_2.txt')
b3 = np.loadtxt('archivostxt/meanaverageV_J10_ro0_Vo-3.txt')
b4 = np.loadtxt('archivostxt/meanaverageV_J10_ro0_Vo-3_2.txt')
# =============================================================================
# Ploteamos las funciones
# Para r
plt.plot(eta1,r)
plt.errorbar(a1[:,0],a1[:,1], a1[:,2], fmt = 'bo')
plt.errorbar(a2[:,0],a2[:,1],a2[:,2],fmt = 'bo')
plt.errorbar(a3[:,0],a3[:,1],a3[:,2],fmt = 'ro')
plt.errorbar(a4[:,0],a4[:,1],a4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','_','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("r")
plt.title("J=10")
plt.xlim([-11,0])
plt.savefig('r_vs_eta_J10.png') 
plt.show()

#Para v
plt.plot(eta1,v)
plt.errorbar(b1[:,0],b1[:,1], b1[:,2], fmt = 'bo')
plt.errorbar(b2[:,0],b2[:,1], b2[:,2], fmt = 'bo')
plt.errorbar(b3[:,0],b3[:,1], b3[:,2], fmt = 'ro')
plt.errorbar(b4[:,0],b4[:,1], b4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','_','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("v")
plt.title("J=10")
plt.xlim([-11,0])
plt.savefig('v_vs_eta_J10.png') 
plt.show()

#%%

# Estudiamos los puntos estables para un J determinado en función de etamedia

J =15
# Primero r en función de eta. 
r = np.linspace(0.04,1.7,100000)

# Sin intensidad de corriente
eta1 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2

# Con intensidad de corriente
eta2 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2 - 3

# Ahora v en función de eta
v = -1/(np.pi*2*r)

# cargamos ficheros de texto con los datos experimentales
a1 = np.loadtxt('archivostxt/meanrate_J15_ro1_Vo0.txt')
a2 = np.loadtxt('archivostxt/meanrate_J15_ro1_Vo0_2.txt')
a3 = np.loadtxt('archivostxt/meanrate_J15_ro0_Vo-3.txt')
a4 = np.loadtxt('archivostxt/meanrate_J15_ro0_Vo-3_2.txt')
# =============================================================================
b1 = np.loadtxt('archivostxt/meanaverageV_J15_ro1_Vo0.txt')
b2 = np.loadtxt('archivostxt/meanaverageV_J15_ro1_Vo0_2.txt')
b3 = np.loadtxt('archivostxt/meanaverageV_J15_ro0_Vo-3.txt')
b4 = np.loadtxt('archivostxt/meanaverageV_J15_ro0_Vo-3_2.txt')
# =============================================================================
# Ploteamos las funciones
# Para r
plt.plot(eta1,r)
plt.errorbar(a1[:,0],a1[:,1], a1[:,2], fmt = 'bo')
plt.errorbar(a2[:,0],a2[:,1],a2[:,2],fmt = 'bo')
plt.errorbar(a3[:,0],a3[:,1],a3[:,2],fmt = 'ro')
plt.errorbar(a4[:,0],a4[:,1],a4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','_','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("r")
plt.title("J=15")
plt.xlim([-11,0])
plt.savefig('r_vs_eta_J15.png') 
plt.show()

#Para v
plt.plot(eta1,v)
plt.errorbar(b1[:,0],b1[:,1], b1[:,2], fmt = 'bo')
plt.errorbar(b2[:,0],b2[:,1], b2[:,2], fmt = 'bo')
plt.errorbar(b3[:,0],b3[:,1], b3[:,2], fmt = 'ro')
plt.errorbar(b4[:,0],b4[:,1], b4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','_','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("v")
plt.title("J=15")
plt.xlim([-11,0])
plt.savefig('v_vs_eta_J15.png') 
plt.show()

#%%

# Estudiamos los puntos estables para un J determinado en función de etamedia

J =20
# Primero r en función de eta. 
r = np.linspace(0.04,2.02,100000)

# Sin intensidad de corriente
eta1 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2

# Con intensidad de corriente
eta2 = -1/(np.pi**2*r**2*4)-J*r+np.pi**2*r**2 - 3

# Ahora v en función de eta
v = -1/(np.pi*2*r)

# cargamos ficheros de texto con los datos experimentales
a1 = np.loadtxt('archivostxt/meanrate_J20_ro1_Vo0.txt')
a2 = np.loadtxt('archivostxt/meanrate_J20_ro1_Vo0_2.txt')
a3 = np.loadtxt('archivostxt/meanrate_J20_ro0_Vo-3.txt')
a4 = np.loadtxt('archivostxt/meanrate_J20_ro0_Vo-3_2.txt')
# =============================================================================
b1 = np.loadtxt('archivostxt/meanaverageV_J20_ro1_Vo0.txt')
b2 = np.loadtxt('archivostxt/meanaverageV_J20_ro1_Vo0_2.txt')
b3 = np.loadtxt('archivostxt/meanaverageV_J20_ro0_Vo-3.txt')
b4 = np.loadtxt('archivostxt/meanaverageV_J20_ro0_Vo-3_2.txt')
# =============================================================================
# Ploteamos las funciones
# Para r
plt.plot(eta1,r)
plt.errorbar(a1[:,0],a1[:,1], a1[:,2], fmt = 'bo')
plt.errorbar(a2[:,0],a2[:,1],a2[:,2],fmt = 'bo')
plt.errorbar(a3[:,0],a3[:,1],a3[:,2],fmt = 'ro')
plt.errorbar(a4[:,0],a4[:,1],a4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','_','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("r")
plt.title("J=20")
plt.xlim([-11,0])
plt.savefig('r_vs_eta_J20.png') 
plt.show()

#Para v
plt.plot(eta1,v)
plt.errorbar(b1[:,0],b1[:,1], b1[:,2], fmt = 'bo')
plt.errorbar(b2[:,0],b2[:,1], b2[:,2], fmt = 'bo')
plt.errorbar(b3[:,0],b3[:,1], b3[:,2], fmt = 'ro')
plt.errorbar(b4[:,0],b4[:,1], b4[:,2], fmt = 'ro')

plt.legend(('teorico','ro=1,Vo=0','_','ro=0,Vo=-3'))
plt.xlabel("eta")
plt.ylabel("v")
plt.title("J=20")
plt.xlim([-11,0])
plt.savefig('v_vs_eta_J20.png') 
plt.show()


# =============================================================================
# #%%
# =============================================================================
# #%%
# # Realizamos el campo vectorial junto con las nullclinas
# etamedia = -5
# J = 15
# I = 0
# 
# # Campo vectorial
# v,r = np.meshgrid(np.linspace(-2.2,1.5,20),np.linspace(0.05,2,20))
# rprima = 1/np.pi + 2*r*v
# vprima = v**2+etamedia+J*r-np.pi**2*r**2+I
# plt.quiver(r,v,rprima,vprima)
# 
# # Nulclina de r
# v = np.linspace(-2.2,1.5,200)
# r = -1/(np.pi*v*2)
# plt.plot(r,v)
# 
# # Nulclina de v
# r = np.linspace(-0.05,2,200)
# v = -np.sqrt(-etamedia - J*r + np.pi**2*r**2)
# plt.plot(r,v)
# plt.xlim([0.05,2])
# plt.ylim([-2.2,1.5])
# plt.xlabel('r')
# plt.ylabel('v')
# plt.legend(('r-nul','v-nul',' '))
# #plt.savefig('campovectorial_J15_eta-5.png') 
# plt.show()
# =============================================================================
