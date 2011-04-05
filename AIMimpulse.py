import numpy as np
import csv
import matplotlib.pyplot as plt

#impulse=np.loadtxt("G:\MathCode\\rssmodelA.txt", delimiter='\t',skiprows=0)
impulseA=np.loadtxt("G:\MathCode\\rssmodelA.txt",delimiter='\t',skiprows=0)
impulseG=np.loadtxt("G:\MathCode\\rssmodelG.txt",delimiter='\t',skiprows=0)
impulseI=np.loadtxt("G:\MathCode\\rssmodelI.txt",delimiter='\t',skiprows=0)
impulseAAA=np.loadtxt("G:\MathCode\\rssmodelAAA.txt",delimiter='\t',skiprows=0)

#impulseA=np.loadtxt("/home/bart/Documents/MathCode/rssmodelA.txt",delimiter='\t',skiprows=0)
#impulseG=np.loadtxt("/home/bart/Documents/MathCode/rssmodelG.txt",delimiter='\t',skiprows=0)
#impulseI=np.loadtxt("/home/bart/Documents/MathCode/rssmodelI.txt",delimiter='\t',skiprows=0)
#impulseA=np.loadtxt("C:\\Users\\bbaker\\Documents\\MacroProg\\Data\\rssmodelA.txt",delimiter='\t',skiprows=0)
#impulseG=np.loadtxt("C:\\Users\\bbaker\\Documents\\MacroProg\\Data\\rssmodelG.txt",delimiter='\t',skiprows=0)
#impulseI=np.loadtxt("C:\\Users\\bbaker\\Documents\\MacroProg\\Data\\rssmodelI.txt",delimiter='\t',skiprows=0)

reader=csv.reader(open("G:\MathCode\impulsetestnames.csv",'rb'))
#reader=csv.reader(open("/home/bart/Documents/MathCode/impulsetestnames.csv",'rb'))
#reader=csv.reader(open("C:\\Users\\bbaker\\Documents\\MacroProg\\Data\\impulsetestnames.csv",'rb'))


names=[]
for x in reader:
    names.append(x[0])
    
sssubtA=np.ones([np.size(impulseA,0),np.size(impulseA,1)])
sssubtG=np.ones([np.size(impulseG,0),np.size(impulseG,1)])
sssubtI=np.ones([np.size(impulseI,0),np.size(impulseI,1)])
sssubtAAA=np.ones([np.size(impulseAAA,0),np.size(impulseAAA,1)])


ss=np.loadtxt("G:\MathCode\\steadystate.txt",delimiter='\t',skiprows=0)
#ss=np.loadtxt("/home/bart/Documents/MathCode/steadystate.txt",delimiter='\t',skiprows=0)
#ss=np.loadtxt("C:\\Users\\bbaker\\Documents\\MacroProg\\Data\\steadystate.txt",delimiter='\t',skiprows=0)

i=0
while i < np.size(sssubtA,1):
    sssubtA[:,i][None]=ss[None]
    i+=1
i=0
while i < np.size(sssubtG,1):
    sssubtG[:,i][None]=ss[None]
    i+=1    
i=0
while i < np.size(sssubtI,1):
    sssubtI[:,i][None]=ss[None]
    i+=1    
i=0
while i < np.size(sssubtAAA,1):
    sssubtAAA[:,i][None]=ss[None]
    i+=1    

plotterA=impulseA-sssubtA
plotterG=impulseG-sssubtG
plotterI=impulseI-sssubtI
plotterAAA=impulseAAA-sssubtAAA


#shock occurs in period 50
#how many fututre periods to simulate forward
futpers=35

plotterA=plotterA[:,50:(49+futpers)]
plotterG=plotterG[:,50:(49+futpers)]
plotterI=plotterI[:,50:(49+futpers)]
plotterAAA=plotterAAA[:,50:(49+futpers)]


#create tech-shock graphs for output and termprem

shty=['A','G','I','AAA']
for m, shock in enumerate(shty):
    i=0
    lengthvec=np.zeros([1,vars()["plotter"+shock].shape[1]])
    zeroguy=np.zeros([1,vars()["plotter"+shock].shape[1]])
    while i < vars()["plotter"+shock].shape[1]:
        lengthvec[0,i]=i+1
        i=i+1
    fig = plt.figure(m)
    if shock !="AAA":
        figtit=plt.suptitle("Shock to "+shock)
    if shock =="AAA":
        figtit=plt.suptitle("Expansionary period: 10 quarters")
    adj = plt.subplots_adjust(hspace=0.6,wspace=0.4)
    j=0
    while j < len(names):
        if names[j]=='Y':
            Yindex=j
        if names[j]=="termprem":
            termpremindex=j
        if names[j]=="Intr":
            intindex=j
        j=j+1
    
    figsub = plt.subplot(311)        
    figsubtit1 = plt.title('Output')
    aaa = plt.plot(lengthvec[0,:],vars()["plotter"+shock][Yindex,:],'b',lengthvec[0,:],zeroguy[0,:],'k')
    figsub = plt.subplot(312)
    figsubtit1 = plt.title('term premium')
    aaa = plt.plot(lengthvec[0,:],vars()["plotter"+shock][termpremindex,:],'b',lengthvec[0,:],zeroguy[0,:],'k')
    figsub = plt.subplot(313)
    figsubtitl = plt.title('real i')
    aaa = plt.plot(lengthvec[0,:],vars()["plotter"+shock][intindex,:],'b',lengthvec[0,:],zeroguy[0,:],'k')
    plt.savefig("G:\\Econ 633\\rssshock"+shock+".png")    
    plt.close
    #plt.savefig("/home/bart/Documents/MathCode/rssshock"+shock+".png")
    #plt.savefig("C:\\Users\\bbaker\\Documents\\MacroProg\\Data\\rssshock"+shock+".png")
    


















