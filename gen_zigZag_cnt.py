# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 16:03:29 2018

@author: OEM
"""

import numpy as np

def zzCNT_gen(n,t,u):
    # see .pdf
    #n chirality
    n = int(n)
    H11 = np.zeros([2*n,2*n]) #Hamiltonian for the supercell
    H12 = np.zeros([2*n,2*n])#coupling to next supercell
    H1 = np.zeros([n,n])# coupling in one A(B) type stack is diagonal in NN-approx
    t1 = np.zeros([n,n])# coupling between A(B) and next B(A) stack on the right (one coupling)
    t2 = np.zeros([n,n])# coupling between B(A) and A(B) 
    #########
    t2[0,n-1] = t #periodic boundary on circumference (closed tube)
    for i in range(n):
        H1[i,i] = u
        t1[i,i] = t
        t2[i,i] = t
    #elements on lower sub-diagonal in t2
    for i in range(n-1):
        t2[i+1,i] = t
     
    #fill supercell H11
    H11[0:n,0:n] = H1
    H11[n:2*n,n:2*n] = H1
    H11[0:n,n:2*n] = t2
    H11[n:2*n,0:n] = np.matrix.transpose(t2)
    #coupling to next supercell        
    H12[n:2*n,0:n] = t1
    
    return H11,H12
    
    
    
    
