# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:32:35 2010

This file is a collection of outside utils

@author: Bart Baker
"""

import numpy as np
import scipy.linalg as linalg

"""The following code is from Alan Isaac's econpy code"""

def hpfilter(y, penalty=1600):
	"""Return: (t,d) trend and deviation
	Based on Lubuele's GAUSS code:
	http://www.american.edu/academic.depts/cas/econ/gaussres/timeseri/hodrick.src
	which is a translation of Prescott's fortran code:
	http://dge.repec.org/codes/prescott/hpfilter.for (Prescott's code)
	which is a 2-pass Kalman filter.
	assumes y is 1d
	penalty is often called 'lambda'

	Conceptualize the calculation as follows:
	eye <- diag( length(y) )
	d2 <- diff( eye, d=2 )
	z <- solve( eye + penalty * crossprod(d2),  y )
	"""
	s = penalty
	assert (s>0)
	n = len(y)
	assert (n>3)
	t = [0]*n  #1d
	d = [0]*n  #1d
	v = np.zeros( (n,3) )
	m1 = y[1]   #changed to zero-based indexing
	m2 = y[0]   #changed to zero-based indexing
	i1 = 3
	i2 = n

	#initialize v
	v11 = 1.0
	v22 = 1.0
	v12 = 0.0
	i = i1    #i initially 3, increments each pass
	istep = 1
	while (i <= i2): #first pass
		#subroutine_pass()
		x = m1
		m1 *= 2
		m1 -= m2
		m2 = x
		x = v11
		z = v12
		v11 = 1/s+4*(x-z)+v22
		v12 = 2*x-z
		v22 = x
		dett = v11*v22-v12*v12
		if istep == 1:  #counting fwd
			v[i-2,0] = v22/dett #changed to zero-based indexing
			v[i-2,2] = v11/dett #changed to zero-based indexing
			v[i-2,1] = -v12/dett #changed to zero-based indexing
			t[i-2] = v[i-2,0]*m1+v[i-2,1]*m2 #changed to zero-based indexing
			d[i-2] = v[i-2,1]*m1+v[i-2,2]*m2 #changed to zero-based indexing
		elif i >= 2: #counting backward
			b11 = v11/dett
			b12 = -v12/dett
			b22 = v22/dett
			e1 = b11*m2+b12*m1+t[i-1] #changed to zero-based indexing
			e2 = b12*m2+b22*m1+d[i-1] #changed to zero-based indexing
			b12 += v[i-1,1] #changed to zero-based indexing
			b22 += v[i-1,2] #changed to zero-based indexing
			b11 += v[i-1,0] #changed to zero-based indexing
			dett = b11*b22-b12*b12
			t[i-1] = (-b12*e1+b11*e2)/dett
		x = v11+1
		z = (y[i-1]-m1)/x #changed to zero-based indexing
		m1 += v11*z
		m2 += v12*z
		z = v11
		v11 -= v11*v11/x
		v22 -= v12*v12/x
		v12 -= z*v12/x
		i += istep
	t[-1] = m1  #ok
	t[-2] = m2  #ok
	m1 = y[-2]  #ok
	m2 = y[-1]  #ok
	i1 = n-2
	i2 = 1
	v11 = 1.0
	v22 = 1.0
	v12 = 0.0
	i = i1
	istep = -1
	while (i >= i2): #second backward pass
		#subroutine_pass()
		x = m1
		m1 *= 2
		m1 -= m2
		m2 = x
		x = v11
		z = v12
		v11 = 1/s+4*(x-z)+v22
		v12 = 2*x-z
		v22 = x
		dett = v11*v22-v12*v12
		if istep == 1:  #counting fwd
			v[i-2,0] = v22/dett #changed to zero-based indexing
			v[i-2,2] = v11/dett #changed to zero-based indexing
			v[i-2,1] = -v12/dett #changed to zero-based indexing
			t[i-2] = v[i-2,0]*m1+v[i-2,1]*m2 #changed to zero-based indexing
			d[i-2] = v[i-2,1]*m1+v[i-2,2]*m2 #changed to zero-based indexing
		elif i >= 2: #counting backward
			b11 = v11/dett
			b12 = -v12/dett
			b22 = v22/dett
			e1 = b11*m2+b12*m1+t[i-1] #changed to zero-based indexing
			e2 = b12*m2+b22*m1+d[i-1] #changed to zero-based indexing
			b12 += v[i-1,1] #changed to zero-based indexing
			b22 += v[i-1,2] #changed to zero-based indexing
			b11 += v[i-1,0] #changed to zero-based indexing
			dett = b11*b22-b12*b12
			t[i-1] = (-b12*e1+b11*e2)/dett
		x = v11+1
		z = (y[i-1]-m1)/x #changed to zero-based indexing
		m1 += v11*z
		m2 += v12*z
		z = v11
		v11 -= v11*v11/x
		v22 -= v12*v12/x
		v12 -= z*v12/x
		i += istep
	t[0] = m1   #changed to zero-based indexing
	t[1] = m2   #changed to zero-based indexing
	i = 0   #changed to zero-based indexing
	while i < n:   #changed to zero-based indexing
		d[i] = y[i]-t[i]
		i = i+1
	return (t,d)

def myHPfilter(y,w):
    t=len(y)
    a=6*w+1
    b=-4*w
    c=w
    d=np.array([[c,b,a]])
    d=np.dot(np.ones([t,1]),d)
    m=np.diag(d[:,2])+np.diag(d[0:t-1,1],1)+np.diag(d[:t-1,1],-1)
    m=m+np.diag(d[:t-2,0],2)+np.diag(d[:t-2,0],-2)
    m[0,0]=1+w
    m[0,1]=-2*w
    m[1,0]=-2*w
    m[1,1]=5*w+1
    m[t-2,t-2]=5*w+1
    m[t-2,t-1]=-2*w
    m[t-1,t-2]=-2*w
    m[t-1,t-1]=1+w
    s=np.dot(linalg.inv(m),y)
    return s
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



