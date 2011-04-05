# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 10:06:58 2010

@author: Barton Baker
"""

#This code is meant to replicate the discretization DSGE solution method found in Alogrithm 2.2 of Fabio Canova's "Methods for Applied Macroeconomics"

import numpy as np
import scipy.linalg as linalg
import discMarkov



def maxinarray(array):
    rows=np.size(array,0)
    cols=np.size(array,1)
    i=0
    max=-5000000000000
    maxrowindex=rows+1
    maxcolindex=cols+1
    while i<rows:
        j=0
        while j<cols:
            if array[i,j]>max:
                max=array[i,j]
                maxrowindex=i
                maxcolindex=j
            j+=1
        i+=1
    if max==-5000000000000:
        return "Error"
    return maxrowindex,maxcolindex,max
                
            



#initialize v matrix

def rbcchoice(endog,exog,origend,origex,V,p):
    #from this we should get all possible utility combinations from in the discrete grid
    #endog should be a vector of possible endognous values
    #exog should be a vector of exogous values
    #origend is the starting point of the endogenous variable
    #origex is the starting point of the exognous variable
    
    #parameters
    ty=0.1
    delta=0.1
    beta=0.9
    phi=2.0
    eta=0.67
    
    outputlist=[]
    i=0
    while i<len(endog):
        #need to create sums of expected utilities times probabilities
        fut=0
        j=0
        while j<np.size(V,1):
            fut=fut+p[origex,j]*V[i,j]
            j+=1
        outputlist=np.append(outputlist,(((1-ty)*endog[origend]**(1-eta)+(1-delta)*endog[origend]-endog[i]-exog[origex])**(1-phi))/(1-phi)+beta*(fut))
        i+=1
    return outputlist



#K=[4.9,5.0,5.075,5.1,5.125,5.2,5.3,5.5,5.7,5.9,6.1,6.2,6.4]
#K=[5.25,5.3,5.35,6.4]
K=[5.3,6.4]
#K=[5.3,5.8,6.4]
#K=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]


#Instead of defining explicit values for the government shock, let us use a user defined method to generate the possible values and Markov transition matrix
gARout=discMarkov.gaussAR(5,1,0.9,1)

G=gARout[0]
#G=[1.1,0.9]

#Also get transition probabilities matrix from gaussAR
#p=np.array([[0.8,0.2],[0.3,0.7]])
p=gARout[1]


#Initialize V matrix as matrix of zeros
V=np.zeros([len(K),len(G)])
#V=np.ones([len(K),len(G)])*-100000

iota=0.01
counter=0
while True:
    #print V
    Vnew=np.zeros([len(K),len(G)])
    i=0
    while i < np.size(V,0):
        j=0
        while j < np.size(V,1):
            Vnew[i,j]=np.max(rbcchoice(K,G,i,j,V,p))
            j+=1
        i+=1
    if np.max(np.abs(Vnew-V)) < iota:
        #print "The limit is the following matrix: "
        #print Vnew
        print "The best estimate of the solution for K and G (with " +str(len(K))+" K grid points and "+ str(len(G))+" G grid points) is:"
        print "K= "+str(K[maxinarray(Vnew)[0]])
        print "G= " +str(G[maxinarray(Vnew)[1]])
        break
    if counter >=50000:
        print "Did not converge"
        break
    counter+=1
    V=Vnew
    
    
#here we will investigate cubic spline interpolation

