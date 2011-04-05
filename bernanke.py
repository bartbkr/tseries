# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:29:40 2010

@author: Bart Baker
"""

import numpy as np
#import utils as util
#import CEEplot
import identVAR

#data=np.loadtxt("G:\MacroProg\Data\VARbernankedata.txt", delimiter='\t', skiprows=1)

#data=np.loadtxt("C:\Users\\bbaker\Documents\MacroProg\Data\VARbernankedata.txt", delimiter='\t', skiprows=1)

#filterdata=data[:,1]
#(hpout,diff)=util.hpfilter(filterdata,14400)
#myfiltdata=util.myHPfilter(filterdata,14400)

#CEEplot.simplot(myfiltdata)
#np.savetxt("C:\Users\\bbaker\Documents\MacroProg\Data\myfiltdata.txt",myfiltdata,delimiter='\t')




"""
datanoseas=np.loadtxt("G:\MacroProg\Data\emplnotseas.txt", delimiter='\t', skiprows=0)
replacedata=[0]*len(filterdata)
i=0
j=0
while i <len(datanoseas):
    if ((i+1) % 13!=0 ):
        replacedata[j]=datanoseas[i]
        i+=1
        j+=1
    if (i!=0 and (i+1) % 13==0):
        i+=1
np.savetxt("G:\MacroProg\Data\cleanedempl.txt",replacedata,delimiter='\t')

(hpout1,diff1)=util.hpfilter(replacedata,14400)
CEEplot.simplot(hpout1)
"""

impdata=np.loadtxt("G:\MacroProg\Data\impberndata.txt",delimiter='\t', skiprows=1)

#FLIP DATA!!!

impdata=np.flipud(impdata)
#get into bernanke form

dates=impdata[:,0]

impdata=impdata[:,1:]

#newdata=np.zeros([np.size(impdata,0)-3,np.size(impdata,1)*4])

lags=4
#i=0
#while i <lags:
#    if i <lags-1:
#        newdata[:,i*np.size(impdata,1):(i+1)*np.size(impdata,1)]=impdata[i:-lags+i+1,:]
#    if i==lags-1:
#        newdata[:,i*np.size(impdata,1):(i+1)*np.size(impdata,1)]=impdata[i:,:]
#    i+=1

#np.savetxt("G:\MacroProg\Data\\tester.txt",newdata,delimiter='\t')

f=open("G:\MacroProg\Data\impberndata.txt",'r')
names=f.readline()
names=names.split()

logs=[4]

VARout=identVAR.VAR(impdata,logs,'Y',lags,'N',25,3,'N','N','N')

betas=VARout['betas']

#get betas into same form as bernanke

prephi=betas.T

constants=prephi[:,0,None]

constantsbank=np.append(constants,np.zeros([(lags-1)*np.size(impdata,1),1]),0)

phi=prephi[:,1:]

addermat=np.eye((lags-1)*np.size(impdata,1))

addermat=np.append(addermat,np.zeros([(lags-1)*np.size(impdata,1),np.size(impdata,1)]),1)

phibank=np.append(phi,addermat,0)

#np.savetxt("G:\MacroProg\Data\\tester.txt",phi,delimiter='\t')


#also get residuals

resids=VARout['resids']



useresids=np.append(resids,np.zeros([(lags-1)*np.size(impdata,1),np.size(resids,1)]),0)

sigmabank=np.dot(useresids,useresids.T)

sigma=np.dot(resids,resids.T)

delta0bank=np.zeros([lags*np.size(impdata,1),1],None)
delta1bank=np.zeros([lags*np.size(impdata,1),1],None)
delta1bank[3,0]=1

delta0=np.zeros([np.size(impdata,1),1],None)
delta1=np.zeros([np.size(impdata,1),1],None)
delta1[3,0]=1


np.savetxt("G:\MacroProg\Data\\bernconstants.txt",constants,delimiter='\t')
np.savetxt("G:\MacroProg\Data\\bernsigma.txt",sigma,delimiter='\t')
np.savetxt("G:\MacroProg\Data\\bernphi.txt",phi,delimiter='\t')
np.savetxt("G:\MacroProg\Data\\berndelta0.txt",delta0,delimiter='\t')
np.savetxt("G:\MacroProg\Data\\berndelta1.txt",delta1,delimiter='\t')













