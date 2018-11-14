# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 14:54:38 2018

@author: OEM
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import MO

def D(n,E,eta,eps,t,sigma):
    h = np.zeros([n,n],dtype=complex)
    h[0,0] = sigma + eps
    h[n-1,n-1] = sigma + eps
    if n > 1:  
        h[n-1,n-2] = t
        h[n-2,n-1] = t
        h[0,1] =     t
        h[1,0] =     t
        for i in range(1,n-1):
            h[i,i] = eps
            if i < n-1:
                h[i,i+1] = t
                h[i+1,i] = t
    return ((E+1j*eta)*np.eye(n)-h) 

def Dblock(n,b,E,eta,eps,t,sigma):
    h = np.zeros([b*n,b*n],dtype=complex)
    h[0:b,0:b] = sigma + eps
    h[(n-1)*b:n*b,(n-1)*b:n*b] = sigma + eps
    if n > 1:  
        h[(n-1)*b:(n)*b,(n-2)*b:(n-1)*b] = t
        h[(n-2)*b:(n-1)*b,(n-1)*b:(n)*b] = t
        h[0:b,b:2*b] =     t
        h[b:2*b,0:b] =     t
        for i in range(1,n-1):
            h[i*b:(i+1)*b,i*b:(i+1)*b] = eps
            if i < n-1:
                h[i*b:(i+1)*b,(i+1)*b:(i+2)*b] = t
                h[(i+1)*b:(i+2)*b,i*b:(i+1)*b] = t
    return  ((E+1j*eta)*np.eye(n*b)-h)  #((E+1j*eta)*np.eye(n*b)-h) 


N = 5
b = 1    
#constructing Hamiltonian for periodic system
#################################################
h11 = 0*np.ones(b)
#s11 = np.ones(b)
#s12 = np.zeros(b)
h12 = 0.5*np.ones(b)
H11 = np.zeros([b,b*N],dtype=np.complex)
S11 = np.ones([b,b*N],dtype=np.complex)
H12 = np.zeros([b,b*(N-1)],dtype=np.complex)
S12 = np.zeros([b,b*(N-1)],dtype=np.complex)
for i in range(N):
    H11[0:b,i*b:(i+1)*b] = h11
    H12[0:b,i*b:(i+1)*b] = h12  
eta = np.array([1e-3],dtype=np.complex)
eps = np.array([1e-3],dtype=np.complex)
E = np.arange(-5,5,0.1,dtype=np.complex)
L = len(E)
dos = np.zeros(E.shape,dtype=np.complex)
sigma = np.zeros(E.shape,dtype=np.complex)
#################################################
for i in range(L):
    sigma[i] = MO.sancho_scalar(E[i],eta,h11,h12,1,0,eps)
    H11[0:b,0:b] = h11 + sigma[i]
    H11[0:b,(N-1)*b:N*b] = h11 + sigma[i]
    g = MO.RGFM(E[i],eta,H11,H12,S11,S12)
    dos[i] = -np.imag(np.sum(g))

plt.plot(E,dos)
plt.xlabel('E')
plt.ylabel('dos')
plt.title('N= '+str(N))
plt.savefig('C:/Users/OEM/Documents/DenizStudium/Masterarbeit/Recherche/Aufgaben/RGFA/Results/DOS_'+str(N)+'.png')
