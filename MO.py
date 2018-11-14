# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:04:40 2018

@author: OEM

Matrix-operations
"""

import numpy as np

#simple recursion expression for periodic systems
##########################################################################################
def gs_rec(E,eta,H_00,H_01,niterations):
#    assert np.imag(E) != 0
    try:  
        Gs = np.linalg.inv((E+1j*eta)*np.eye(len(H_00))-H_00) 
        for i in range(niterations):
            Gs = np.linalg.inv(E*np.eye(len(H_00))-H_00-np.matrix.getH(H_01)@Gs@H_01)
    except TypeError:
        Gs = 1/((E+1.j*eta)-H_00)
        for i in range(niterations):
            Gs = 1./(E+1.j*eta-H_00-np.conjugate(H_01)*Gs*H_01)
        sigma = np.conjugate(H_01)*Gs*H_01
    return sigma  
##########################################################################################       

#sancho method
##########################################################################################    
def sancho(energy,h,t0_matrix,sh,st,eps):
    es = energy*sh-h
    e = energy*sh-h
    a = energy*st-t0_matrix
    b = energy*np.matrix.getH(st) - np.matrix.getH(t0_matrix)
    

    
    
    
    while((np.linalg.norm(abs(a), ord='fro')) > eps):
        g = np.linalg.inv(e)
        bga = b @ g @ a
        agb = a @ g @ b
        e = e - bga - agb
        es = es - agb

        a = -a @ g @ a
        b = -b @ g @ b

    G = np.linalg.inv(es)
    SIGMA = np.matrix.getH(t0_matrix)@G@t0_matrix
    return SIGMA
##########################################################################################
def sancho_scalar(energy,eta,h,t,sh,st,eps):
    es =np.array([(energy+1j*eta)*sh-h],dtype=np.complex)
    e = np.array([(energy+1j*eta)*sh-h],dtype=np.complex)
    a = np.array([(energy+1j*eta)*st-t],dtype=np.complex)
    b = np.array([(energy+1j*eta)*np.conjugate(st)-np.conjugate(t)],dtype=np.complex)
    while(abs(a) > eps):
        g = 1./e
        bga = b*g*a
        agb = a*g*b
        e = e - bga -agb
        es = es - agb
        
        a = -a*g*a
        b = -b*g*b
    sigma = 1./es*t**2

    return sigma  
##########################################################################################
#Recursive GF Method
def rgfm(E,eta,N,H11,H12,S11,S12,Sigma):
    #periodic system only diagonal elements
    #H11, H12, S11, S12 are bxb matrices
    H11 = np.asmatrix(H11,dtype=np.complex)
    H12 = np.asmatrix(H12,dtype=np.complex)
    S11 = np.asmatrix(S11,dtype=np.complex)
    S12 = np.asmatrix(S12,dtype=np.complex)
    E = np.array([E],dtype=np.complex)
    b = H11.shape[0]
    g = np.zeros([b*N,b],dtype=np.complex)
    G = np.zeros([b*N,b],dtype=np.complex)
    g[0:b,0:b] = np.linalg.inv((E+1j*eta)*np.ones(b)-H11-Sigma)
#    G[(N-1)*b:N*b,0:b] = np.linalg.inv((E+1j*eta)*np.ones(b)-H11-Sigma)
#    

    for i in range(1,N):
        g[i*b:(i+1)*b,0:b]=np.linalg.inv((E+1j*eta)*np.ones(b,dtype=np.complex)-H11-np.matrix.getH(H12)@g[(i-1)*b:i*b,0:b]@H12)  #g OK!!!!
    
    G[(N-1)*b:N*b,0:b] = np.linalg.inv((E+1j*eta)*np.ones(b,dtype=np.complex)-H11-Sigma)         #g[(N-1)*b:N*b,0:b]
    
    for i in range(1,N):
        G[(N-1-i)*b:(N-i)*b,0:b] = g[(N-1-i)*b:(N-i)*b,0:b]@(np.eye(b,dtype=np.complex)+H12@G[(N-i)*b:(N-i+1)*b,0:b]@np.matrix.getH(H12)@g[(N-1-i)*b:(N-i)*b,0:b])
    return G


def RGFM(E,eta,H11,H12,S11,S12):
    #Self-Energy has to be inseted in input.
    #INPUT: H11 = [H11[0:b,0:b],...,H[0:b,(k-1)*b:k*b],...,H[0:b,(N-1)*b:N*b]] 0th element, kth element N-1th element
    b = H11.shape[0] #blocksize
    N = int(H11.shape[1]/b) #size of central cystem
    A11 = (E+1j*eta)*S11-H11
    A12 = (E+1j*eta)*S12-H12#the other 'pseudo-diagonal' is the same with all elements transposed and conj.
    gdiagL = np.zeros([b,N*b],dtype=np.complex)#init diagonal g=[g[0:b,0:b],...,g[0:b,(k-1)*b:k*b],...,g[0:b,(N-1)*b:N*b]]
    gdiagR = np.zeros([b,N*b],dtype=np.complex)
    Gdiag = np.zeros([b,N*b],dtype=np.complex)#init diagonal G same as g
    Gi1 = np.zeros([b,N*b],dtype=np.complex)
    GiN = np.zeros([b,N*b],dtype=np.complex)
    #forward left#################################
    gdiagL[0:b,0:b] = np.linalg.inv(A11[0:b,0:b])
    for i in range(1,N):
        gdiagL[0:b,i*b:(i+1)*b]=np.linalg.inv(A11[0:b,i*b:(i+1)*b]-np.matrix.getH(A12[0:b,(i-1)*b:i*b])@gdiagL[0:b,(i-1)*b:i*b]@A12[0:b,(i-1)*b:i*b]) 
    
    #backward left#################################
    #last element of Gdiag is equal to last element of gdiag
    Gdiag[0:b,(N-1)*b:N*b] = gdiagL[0:b,(N-1)*b:N*b]
    
    for i in range(1,N):
        Gdiag[0:b,(N-1-i)*b:(N-i)*b] = gdiagL[0:b,(N-1-i)*b:(N-i)*b]@(np.eye(b)+A12[0:b,(N-1-i)*b:(N-i)*b]@Gdiag[0:b,(N-i)*b:(N+1-i)*b]@np.matrix.getH(A12[0:b,(N-1-i)*b:(N-i)*b])@gdiagL[0:b,(N-1-i)*b:(N-i)*b]) 
     
    #forward R#################################
    gdiagR[0:b,(N-1)*b:N*b] = np.linalg.inv(A11[0:b,(N-1)*b:N*b])
    for i in range(1,N):
        gdiagR[0:b,(N-1-i)*b:(N-i)*b] = np.linalg.inv(A11[0:b,(N-1-i)*b:(N-i)*b]-A12[0:b,(N-1-i)*b:(N-i)*b]@gdiagR[0:b,(N-i)*b:(N+1-i)*b]@np.matrix.getH(A12[0:b,(N-1-i)*b:(N-i)*b]))
    
    #GiN from gdiagL#################################
    GiN[0:b,(N-1)*b:N*b] = Gdiag[0:b,(N-1)*b:N*b]
    for i in range(1,N):
        GiN[0:b,(N-1-i)*b:(N-i)*b] = -gdiagL[0:b,(N-1-i)*b:(N-i)*b]@A12[0:b,(N-1-i)*b:(N-i)*b]@GiN[0:b,(N-1)*b:N*b]
        
        
        
    #Gi1#################################
    Gi1[0:b,0:b] = Gdiag[0:b,0:b]
    for i in range(1,N):
        Gi1[0:b,i*b:(i+1)*b] = -gdiagR[0:b,i*b:(i+1)*b]@A12[0:b,(i-1)*b:(i)*b]@Gi1[0:b,(i-1)*b:(i)*b]
    
    
    
    
    
    return Gdiag

##########################################################################################