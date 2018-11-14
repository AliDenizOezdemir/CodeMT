# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 11:10:45 2018

@author: OEM
"""
import numpy as np
import os
import yaml

def get_matrices(string):
    H = []
    S = []
    f = open(string,'r')
    Y = yaml.load(f.read())
    f.close()
    H_S = np.loadtxt(Y['folder']+'/H_S.txt', usecols=None)     
    H_L = np.loadtxt(Y['folder']+'/H_L.txt', usecols=None)    
    H_R = np.loadtxt(Y['folder']+'/H_R.txt', usecols=None)    
    T_L = np.loadtxt(Y['folder']+'/T_L.txt', usecols=None)   
    T_R = np.loadtxt(Y['folder']+'/T_R.txt', usecols=None)    
    T_LS = np.loadtxt(Y['folder']+'/T_LS.txt', usecols=None)    
    T_RS = np.loadtxt(Y['folder']+'/T_RS.txt', usecols=None)    
    S_S = np.loadtxt(Y['folder']+'/S_S.txt', usecols=None)   
    S_R = np.loadtxt(Y['folder']+'/S_R.txt', usecols=None)   
    S_L = np.loadtxt(Y['folder']+'/S_L.txt', usecols=None)   
    S_LS = np.loadtxt(Y['folder']+'/S_LS.txt', usecols=None)   
    S_RS = np.loadtxt(Y['folder']+'/S_RS.txt', usecols=None)   
    S_SR = np.loadtxt(Y['folder']+'/S_SR.txt', usecols=None)   
    S_SL = np.loadtxt(Y['folder']+'/S_SL.txt', usecols=None)
    #matrices are stored in dictionary.
    return {'H_S': H_S,'H_L': H_L,'H_R': H_R,'T_L': T_L,'T_R': T_R,'T_LS': T_LS,'T_RS': T_RS,'S_S': S_S,'S_R': S_R,'S_L': S_L,'S_LS': S_LS,'S_RS': S_RS,'S_SR': S_SR,'S_SL': S_SL}
    
def write_data(path,filename,file):
    if not os.path.exists(path):
        os.makedirs(path)
    
    np.savetxt(path+'/'+filename + '.txt',file)
    