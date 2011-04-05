# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:47:52 2010

Transition probabilties for Discrete First Order AutoRegressive Process

@author: bbaker
"""
import numpy as np
from scipy.stats import norm

def gaussAR(n,m,lamda,sigmasqr):
    """
    # n is the number of discrete points needed
    # m is multiple of standard deviation to use for extremem values
    # lambda is the autoregressive parameter
    # sigmasqr is the variance of the shock
    
    #returns tuple:
    # (posvals,mark)
    # posvals a list of the discrete possible values
    # mark is the Markov probability transition matrix corresponding to posvals
    
    #Tauchen(1985) method for selecting values 
    """
    #generate possible values by multiplying m by unconditional standard deviation
    max=m*(sigmasqr/(1-lamda**2))**(0.5)
    min=-max
    sigma=np.sqrt(sigmasqr)
    posvals=[min]
    i=1
    while i<n:
        posvals=np.append(posvals,min+((max-min)/(n-1)*i))
        i+=1
    w=posvals[1]-posvals[0]#this will be the step between each point
    #might want to include more points closer to zero
    j=0
    k=0
    mark=np.zeros([n,n])
    while j<n:
        if k==0:
            mark[j,k]=norm.cdf((posvals[0]-lamda*posvals[j]+w/2)/sigma)
            k+=1
        if 1<=k<=n-1:
            mark[j,k]=norm.cdf((posvals[k]-lamda*posvals[j]+w/2)/sigma)-norm.cdf((posvals[k]-lamda*posvals[j]-w/2)/sigma)
            k+=1
        if k==n-1:
            mark[j,k]=1-norm.cdf((posvals[n-1]-lamda*posvals[j]-w/2)/sigma)
            k=0
            j+=1
    return posvals, mark
        