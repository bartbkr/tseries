# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 14:31:24 2010

@author: Barton Baker
"""
import numpy as np
from numpy import *

numvars=9

L_n=np.zeros([(1.0/2)*(numvars)*(numvars+1),numvars**2])#this is the elimination matrix


#this method is taken from Magnus and Neudecker (1980)
ident=np.eye(numvars)
j=0
while j < numvars:
    i=j
    while j <= i < numvars:
        u=np.zeros([(1.0/2)*numvars*(numvars+1),1])#initialize u_ij matrix
        u[(j)*numvars+(i+1)-(1.0/2)*(j+1)*j-1]=1 #because of 0 indexing, had to mess with formula from magnus (1980) paper
        L_n=L_n+kron(kron(u,ident[:,j].T),ident[:,i].T)
        i=i+1
    j=j+1
    
#output['L_n']=L_n #this is the Elimination matrix

#since the variance/covariance matrix is symmetric, the Moore-Penrose inverse of L is D

D_nT=np.zeros([L_n.shape[0],L_n.shape[1]])

j=0
while j < numvars:
    i=j
    while j <= i < numvars:
        u=np.zeros([(1.0/2)*numvars*(numvars+1),1])#initialize u_ij matrix
        u[(j)*numvars+(i+1)-(1.0/2)*(j+1)*j-1]=1 #because of 0 indexing, had to mess with formula from magnus (1980) paper
        Tij=np.zeros([numvars,numvars])
        Tij[i,j]=1
        Tij[j,i]=1
        D_nT=D_nT+dot(u,(Tij.ravel('F')[:,None]).T)
        i=i+1
    j=j+1

D_n=D_nT.T

x=numvars #x keeps track of sum function to take out correct columns
y=1 #y keeps track of lost columns
collector=D_n
collector=np.delete(collector,0,1)
###################here's the problem
while x >= 2:
    collector=np.delete(collector,sum(b for b in range(x,numvars+1))-y,1)
    #print sum(b for b in range(x,numvars+1))
    x=x-1
    y=y+1
    #print collector[:,0]
    #print raw_input("Continue?")
#output['collector']=collector
S_Bpre=collector

#output['S_Bpre']=S_Bpre
        
        
S_B=np.zeros([S_Bpre.shape[0],S_Bpre.shape[1]])
i=0
j=0
while j < S_Bpre.shape[1]:
    while i < S_Bpre.shape[0]:
        if S_Bpre[i,j]==1:
            S_B[i,j]=1
            S_B[i+1:,j][:,None]=np.zeros([S_Bpre.shape[0]-i-1,1])
            break
        i=i+1
    j=j+1     
        
        
S_BT=np.zeros([S_Bpre.shape[0],S_Bpre.shape[1]])
j=0
while j < S_Bpre.shape[1]:
    i=S_Bpre.shape[0]-1
    while i > 0:
        if S_Bpre[i,j]==1:
            S_BT[i,j]=1
            S_BT[:i,j][:,None]=np.zeros([i,1])
            break
        i=i-1
    j=j+1
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
