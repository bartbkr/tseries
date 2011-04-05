import numpy as np
import matplotlib.pyplot as plt

def CEEplotter(impulse,diff,loci,upci,names):
    lengthvec=np.zeros([1,impulse.shape[1]])
    
    i=0
    m=int(np.ceil(np.sqrt(len(names))))*100
    n=int(np.ceil(np.sqrt(len(names))))*10
    
    
    if diff=='Y':
        while i < impulse.shape[1]:
            lengthvec[i,0]=i+1
            i=i+1
        plt.ion()
        fig1 = plt.figure(1)
        adj = plt.subplots_adjust(hspace=0.6,wspace=0.4)
        j=0
        while j < impulse.shape[0]:
            figsub = plt.subplot(m+n+j+1)
            figsubtit1 = plt.title(names[n])
            aaa = plt.plot(lengthvec,np.cumsum(impulse[:,j:j+1]))
            j=j+1
        plt.show
    

    if diff=='N':
        while i < impulse.shape[1]:
    		lengthvec[0,i]=i+1
    		i=i+1
        zeroguy=np.zeros([1,impulse.shape[1]])
    	plt.ion()
    	fig1 = plt.figure(1)
    	adj = plt.subplots_adjust(hspace=0.6,wspace=0.4)
        j=0
        while j < len(names):
            figsub = plt.subplot(m+n+j+1)
            figsubtit1 = plt.title(names[j])
            aaa = plt.plot(lengthvec[0,:],impulse[j,:],'b',lengthvec[0,:],impulse[j,:]+upci[j,:],'c',lengthvec[0,:],impulse[j,:]+loci[j,:],'c',lengthvec[0,:],zeroguy[0,:],'k')
            j=j+1
    	plt.show
    

def CEEplotsimple(impulse,names):
    i=0
    m=int(np.ceil(np.sqrt(len(names))))*100
    n=int(np.ceil(np.sqrt(len(names))))*10
    lengthvec=np.zeros([1,impulse.shape[1]])
    while i < impulse.shape[1]:
        lengthvec[0,i]=i+1
        i=i+1
    zeroguy=np.zeros([1,impulse.shape[1]])
    plt.ion()
    fig1 = plt.figure(1)
    adj = plt.subplots_adjust(hspace=0.6,wspace=0.4)
    j=0
    while j < len(names):
        figsub = plt.subplot(m+n+j+1)
        figsubtit1 = plt.title(names[j])
        aaa = plt.plot(lengthvec[0,:],impulse[j,:],'b',lengthvec[0,:],zeroguy[0,:],'k')
        j=j+1
    plt.show

def simpleplot(impulse,names,var):
    j=0
    while j < len(names):
        if names[j]==var:
            index=j
        j=j+1
    i=0
    lengthvec=np.zeros([1,impulse.shape[1]])
    zeroguy=np.zeros([1,impulse.shape[1]])
    while i < impulse.shape[1]:
        lengthvec[0,i]=i+1
        i=i+1
    plt.ion()
    fig1 = plt.figure(1)
    figtitle=plt.title(var)
    aaa=plt.plot(lengthvec[0,:],impulse[index,:],'b',lengthvec[0,:],zeroguy[0,:],'k')
    plt.show
    
def simplot(impulse):
    """One dimensional"""
    i=0
    lengthvec=np.zeros([1,len(impulse)])
    while i < len(impulse):
        lengthvec[0,i]=i+1
        i=i+1
    plt.ion()
    fig1 = plt.figure(1)
    aaa=plt.plot(lengthvec[0,:],impulse[:])
    plt.show
    
    
    
    
    
    
    
    

