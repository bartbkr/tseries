# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 16:15:08 2010

@author: Lab User
"""
from numpy import *

def vech(A):
    length=A.shape[1]
    a=length
    sum=0
    while 0 < a <= length:
        sum=sum+a
        a=a-1
    
    i=0
    vechvec=[]
    while i < length:
        b=i
        while b < length:
            vechvec.append(A[b,i])
            b=b+1
        i=i+1
    vechvec=array(vechvec)
    return vechvec
