# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 09:35:55 2018

@author: OEM
"""

import MO
import matplotlib.pyplot as plt
import numpy as np
import yaml
import importlib.util
spec = importlib.util.spec_from_file_location("NEGFGlobal", "C:/Users/OEM/Documents/DenizStudium/Masterarbeit/Recherche/Aufgaben/frompramit/NEGFGlobal.py")
foo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(foo)

H,S = foo.ao_file_loader('C:/Users/OEM/Documents/DenizStudium/Masterarbeit/Recherche/Aufgaben/Cp2k/CNT','nt-3-0-1.ao')
b = 13
N = int(H.shape[0]/b)
h11 = H[0:b,0:b]
s11 = S[0:b,0:b]
h12 = H[0:b,b:2*b]
s12 = S[0:b,b:2*b]

h11 = np.array(h11,dtype=np.complex)
s11 = np.array(s11,dtype=np.complex)
s12 = np.array(s12,dtype=np.complex)
h12 = np.array(h12,dtype=np.complex)
H11 = np.zeros([b,N*b],dtype=np.complex)
S11 = np.zeros([b,N*b],dtype=np.complex)
H12 = np.zeros([b,(N)*b],dtype=np.complex)
S12 = np.zeros([b,(N)*b],dtype=np.complex)

for i in range(N):
   
    H11[0:b,i*b:(i+1)*b] = h11
    S11[0:b,i*b:(i+1)*b] = s11
    
for i in range(N-1):
    H12[0:b,i*b:(i+1)*b] = h12
    S12[0:b,i*b:(i+1)*b] = s12

eta = np.array([1e-3],dtype=np.complex)
eps = np.array([1e-3],dtype=np.complex)
E = np.arange(-5,5,0.1,dtype=np.complex)
EE = E+1j*eta
L = len(E)
dos = np.zeros(E.shape,dtype=np.complex) 
sigma = np.zeros([b,b*L],dtype=np.complex)   #sigma[E[i]] = sigma[0:b,i*b:(i+1)*b]
for i in range(L):
    sigma[0:b,i*b:(i+1)*b] = MO.sancho(EE[i],h11,h12,s11,s12,eps) # def sancho(energy,eta,h,t0_matrix,sh,st,eps):
#    H11[0:b,0:b] = h11 + sigma[0:b,i*b:(i+1)*b]
#    H11[0:b,(N-1)*b:N*b] = h11 + sigma[0:b,i*b:(i+1)*b]
#    g = MO.RGFM(E[i],eta,H11,H12,S11,S12)
    dos[i] = -1/np.pi*np.real(np.trace(sigma[0:b,i*b:(i+1)*b]))
    
#change in comment
plt.plot(E,dos)
plt.xlabel('E')
plt.ylabel('Self-energy from Sancho')
plt.savefig('C:/Users/OEM/Documents/DenizStudium/Masterarbeit/Recherche/Aufgaben/Cp2k/CNT/SIGMA.jpg')