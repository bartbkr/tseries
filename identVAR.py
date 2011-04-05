#Barton Baker

#Run VAR with recursive shock identification


import numpy as np
from numpy import *
import scipy as Sci
import scipy.linalg
import CEEplot #custom plotting program
import Pyplus as pyplus

def VAR(data,logs,const,t,diff,macoeffs,shockvar,orthog,showplot,errorbands):
    """
    logs (list): should be a list of the columns of the data file that should be logged before performing analysis (remember 0 indexing)
    const (string): 'Y' or 'N' for adding constant to regression
    t (int): time periods back for VAR (t-n with n periods)
    diff (string): 'Y' or 'N' for first differencing data before performing analysis
    maceoffs (int): number of moving average coefficient matrices to calculate, this is also the number of periods to run the IRF
    shockvar (int): var that will be shocked to test for impact
    orthog (str): 'Y' or 'N' whether to orthognoize (recursive identification) of shocks
    showplot (string): 'Y' or 'N' whether to show plot for IRFs
    errorbands (str): 'Y' or 'N' whether to show error bands
    
    
    created variables:
    names: names of the variables in order
    betas: estimated constants (first column) and autoregressive coefficients which follow in groups of columns by t going further back
    betas2: estimate of betas using alternative method directly copied from Hamilton(1994) page 393 ML
    resids: residuals from ML estimation
    mu: unconditonal mean
    Q: this corresponds to the A matrix from Hamilton(1994) page 320, which is lower triangular and has 1s along the diagonal
    irfn: impulse response function matrix for effect of shock to var n at time t
    psi_s: psi matrix of moving average coefficients for t+s with shock at time t
    A: A matrix as in Hamilton (1994)
    D: D matrix as in Hamilton (1994)
    G_s: iterative matrix required to calculate standard errors using analytical derivatives for t+s
    Q_T: matrix used to cumpute standard errors, constructing by crossing x_t data and dividing by number of obs
    variances_s: these are the variance error matrices for t+s (appear at top of page 337 of Hamilton(1994)), take square root of diagonal element
    L_n: Elimination matrix for calculating Duplication matrix
    D_n: Duplication matrix for calculating standard errors of orthogonalized IRF
    B_0: B naught matrix used in construction of standard errors of orthogonalized IRF
    S_B: S_B matrix from Hamilton for vec(B_0) calculation
    S_BT: S_BT matrix from Hamilton for vec(B_0T) calculation
    S_D: S_D matrix from Hamilton for vec(D) calculation
    J: J matrix from Hamilton
    G_B: G_B matrix from Hamilton
    G_D: G_D matrix from Hamilton
    upci: upper error band for impulse response function
    loci: lower error band for impulse response function
    
    output is stored in dictionary
    """
    output={}
    
    nperiod=macoeffs
    
    
    #nipa1 = np.loadtxt('H:\Data\CEEdata.txt', delimiter='\t', skiprows=1) #import as tab-delimited data file
    nipa1 = data
    #nipa1=np.loadtxt('G:\MacroProg\Data\CEEdata.txt', delimiter='\t', skiprows=1) #import as tab-delimited data file
    
    #we should flip the data matrix to match the literature, so data moves from top to bottom to earlier dates
    #nipa1=np.flipud(nipa)
    
    #we can fetch the variable names
    #f=open("C:\Users\Bart\Documents\Macro Programming\Data\CEEdata.txt",'r')
    #f=open("G:\MacroProg\Data\CEEdata.txt",'r')
    #names=f.readline()
    #names=names.split()
        
    #output['names']=names
    
    
    #these were all specific to the first file (nipa1)
    # log of quarter of data 1965.5-1995.5
    #gdp=nipa1[:,0]# log of real gdp
    #C=nipa1[:,1]# log of real consumption
    #gdpP=nipa1[:,2]# log of gdp deflator
    #I=nipa1[:,3]# log of real investment
    #W=nipa1[:,4]# log of real wage
    #prod=nipa1[:,5]# log of Non-farm business sector productivity
    #irate=nipa1[:,6] # log of interest rate
    #profit=nipa1[:,7]#log real profit
    #m2=nipa1[:,8] #log quaterly rate of change of M2
    
    
    
    #######################This section logs the columns of data that the user requests
    
#    loggers=raw_input("Are there columns that you would like to take the log of? (Y/N):   ")
#    i=0
#    if not (loggers=='Y' or loggers=='N'):
#    	loggers=raw_input("What?:  ")
    
    
#    if loggers=='Y':
    for x in logs:
    	#numnolog=int(raw_input("How many?:   "))
    	#while i<numnolog:
        #x=int(raw_input("Which column does it lie in (remember 0 indexing)?:    "))
        nipa1[:,x]=np.log(nipa1[:,x])
    
    output['nipa1']=nipa1
    
    
    numvars=nipa1.shape[1]
    
    ###########################add constant to VAR???
    
    #const=raw_input("Would you like to add a constant to the VAR?(Y/N):    ")
    
    #if not (const=='Y' or const=='N'):
    #	const=raw_input("What?:  ")
    
    ###########################Queries for t-?
    
    #t=int(raw_input("How many time periods back does the VAR go? (t-n with n periods):    "))
    
    #diff=raw_input("Would you like to first difference the data? (Y/N):    ")
    
    ##########This section gets data into single matrix
    
    #if not (diff=="Y" or diff=="N"):
    #	diff=raw_input("What?:  ")
    
    if diff=="Y":
    	data=np.zeros([len(nipa1)-1,numvars])
    	i=0
    	while i < numvars:# this step just gets data into single matrix without dates, and first differences the data
    		data[:,i]=nipa1[0:-1,i]-nipa1[1:,i]
    		i=i+1
    
    if diff=="N":
    	data=nipa1
    	#i=0
    	#while i < numvars:# this step just gets data into single matrix without dates
    	#	data[:,i]=nipa1[:,i]
    	#	i=i+1
    
    
    #####################now run VAR
    
    #np.savetxt('G:\MacroProg\Data\data.txt',data, delimiter='\t') #just to test output
    
    yt=data[0:-t,] # yts
    yt_=np.zeros([len(yt),t*numvars])
    
    i=0
    while i < t: #this creates xs which are lagged yts by t
        if i< t-1:
            yt_[:,(i*numvars):((i+1)*numvars)]=data[(i+1):((-t+1)+i),:]
            i=i+1
        if i==(t-1):
            yt_[:,(i*numvars):((i+1)*numvars)]=data[(i+1):,:]
            i=i+1
    
    
    #np.savetxt('G:\MacroProg\Output\yt_file.txt', yt_, delimiter='\t')
    
    output['yt']=yt
    output['yt_']=yt_
    
    
    
    
    ####now get betas
    ####we will use two methods to ensure accuracy
    
    ###########betas method 1
    
    
    if const=='Y':
        constmat=np.ones([len(yt_),numvars*t+1])
        constmat[:,1:]=yt_
        betas=dot(linalg.pinv(constmat),yt)
        
    #np.savetxt('G:\MacroProg\Output\constmat.txt',constmat,delimiter='\t')        
        
    if const=='N':
        betas=dot(linalg.pinv(yt_),yt)
    
    
    
    output['betas']=betas
    
    
    
    ###########betas method 2
    
    #this method follows to a t the method performed on page 293 of Hamilton
    
    #let us transpose the yt_matrix and yt matrix to get them in the same form as in hamilton
    
    yt_vert=yt_.T
    
    ytvert=yt.T
    
    y_wcoll=np.zeros([numvars,1])
    
    if const=='Y':
        x_wcoll=np.zeros([numvars*t+1,numvars*t+1]) #create matrix to hold summed values of matrix sums in OLS
        y_wcoll=np.zeros([numvars,numvars*t+1])
        onesie=np.ones([len(yt_vert)+1,yt_vert.shape[1]])#add row of ones to yt_ matrix
        onesie[1:]=yt_vert
        yt_vert=onesie # set first row of matrix to zero
    
    if const=='N':
        x_wcoll=np.zeros([numvars*t,numvars*t])
        y_wcoll=np.zeros([numvars,numvars*t])
    
    
    w=0 # this section sums all of the dot products of the single vectors which includes collections of t periods of data
    
    while w < len(yt_):
        xxt=dot(yt_vert[:,w, None],yt_vert[:,w, None].T) #this is the x's times the x's
        x_wcoll =x_wcoll+xxt
        xyt=dot(ytvert[:,w, None],yt_vert[:,w, None].T) #this is the y's times the x's
        y_wcoll=y_wcoll+xyt
        w=w+1
    
    
    betas2=dot(y_wcoll,linalg.inv(x_wcoll))
    
    output['x_wcoll']=x_wcoll
    
    output['betas2']=betas2
    
    #np.savetxt('G:\MacroProg\Output\mbetas2.txt',betas2,delimiter='\t')
    
    
    ##now onto collecting residuals
    
    #transpose betas matrix for easier irf creation
    
    betasT=betas.T
    
    resids=ytvert-dot(betasT,yt_vert)#collect residuals
    
    output['resids']=resids
    
    #We can use this to compute the implied estimator standard errors
    
    F=np.zeros([numvars*t,numvars*t])
    #F[0:numvars,:]=betas.T
    output['F']=F
    
    
    
    #Unconditional Mean
    i=0
    mupre=np.eye(numvars) #all of the 
    while i < t:
        mupre=mupre-betasT[:,i*numvars+1:(i+1)*numvars+1]
        i=i+1
        
    
    
    mu=dot(linalg.inv(mupre),betasT[:,0, None]) #unconditional mean
    
    output['mu']=mu
    
    

    
    
    
    
    #################At this point, we can use the autoregressive coefficients to estimate the implied moving average coefficients
    
    #macoeffs=int(raw_input("How many MA coefficients would you like to include?:    "))
    
    #nperiod=int(raw_input("How many periods would you like the IRF's to run?:   "))
    
    #nperiod + t should be used so that macoeffs can line up with the irf length
    nperiods=nperiod+t+1
    
      
    #first we create 9 impulse response functions to hold the effects of a one-time shock to one of the variables
    
    d=0
    while d < numvars:
        vars()['irf'+str(d)]=np.zeros([numvars,nperiods])
        d=d+1
        
    #next, we run an impulse response functions, testing the effect of a one-time shock to each of the variable each time
    #notice that we do NOT include the constants in this for we are trying to pull out the implied moving average coefficients
    
    k=0
    while k < numvars:
        vars()['shock'+str(k)]=np.zeros([numvars,nperiods])
        vars()['shock'+str(k)][k,t]=1.0 #change this to t-1?
        i=t
        while i < nperiods:
            c=0
            while c < numvars:
                z=0
                holder=0
                while z < t:
                    holder=holder+dot(betasT[c,numvars*z+1:numvars*(z+1)+1],vars()['irf'+str(k)][:,i-1-z])
                    z=z+1
                vars()['irf'+str(k)][c,i]=holder+vars()['shock'+str(k)][c,i]
                c=c+1
            i=i+1
        k=k+1
    
    #output['holder']=holder
    
    #output['betasT']=betasT
    
      
    #save these to the output set just for checkup
    
    #in this section, we replace the elements of the psi_ matrices with the corresponding values. 
    #The (i,j) element of the psi_s matrix refers to the effect of a one-time shock on the jth variable to the ith variable at time t+s.
    
    #remember psi_1 matrix should refer to the period directly following the shock
    
    s=1
    while s <= macoeffs:
        vars()['psi_'+str(s)]=np.zeros([numvars,numvars])
        j=0
        while j < numvars:
            i=0
            while i < numvars:
                vars()['psi_'+str(s)][i,j]=vars()['irf'+str(j)][i,s+t]
                i=i+1
            j=j+1
        s=s+1
    
    
    #now output IRFs and macoeffs to dictionary
    #IRFs
    
    i=0
    while i < numvars:
        output['irf'+str(i)]=vars()['irf'+str(i)]
        i=i+1
        
    i=0
    while i < numvars:
        output['shock'+str(i)]=vars()['shock'+str(i)]
        i=i+1
        
    
    
  
    
    
    
    ###################################IRF
    
    #we would like to identify the VAR using recursive identification
    
    #create omegahat
    omegahat=dot(resids,resids.T)*(1.0/(len(yt_)-len(betas))) #this corrects for the number of estimated coefficients
    
    output['omegahat']=omegahat
    
    #R matrix (P matrix in Hamilton)
    P=linalg.cholesky(omegahat)
    
    output['P']=P
    
    #np.savetxt('G:\MacroProg\Data\P.txt',P, delimiter='\t')
    
    #the following creates the Q matrix so that phi0 can be computed (Rnondiag is A in Hamilton
    
    Donehalf=P*eye(len(P)) #elementwise multiply eliminates non-diagonal elements
    
    D=dot(Donehalf,Donehalf)
    
    output['D']=D

    A=dot(P,linalg.inv(Donehalf))
    
    output['A']=A #this is the equivalent of the A matrix from Hamilton
    
    
    ##We can use this information to calculate the standard errors for the IRFs
    ##notation corresponds to Hamilton
    
    
    pi=betas.ravel('F')[:,None] #this is a single vector of the constants and autoregressive coefficients
    
    output['pi']=pi
    
    i=1
    while i <= macoeffs:
        vars()['lilpsi_'+str(i)]=np.ravel(vars()['psi_'+str(i)].T, order='F')
        i=i+1
        
    s=1    
    while s <=macoeffs:
        output['lilpsi_'+str(s)]=vars()['lilpsi_'+str(s)]
        s=s+1
    
    
    
    
    
    ###########the first method we can use to determine standard errors is based on analytical derivatives (Hamilton(1994) page 336)
    ####this method is only appropriate for non-orthogonized IRFs
    
    #set all psi_(-) to arrays of zeros and psi_0 to identity matrix
    
    i=1
    while i <=t:
        vars()['psi_'+str(-i)]=np.zeros([numvars,numvars])
        i=i+1
        
    psi_0=eye(numvars)
    
    
    s=-t
    while s <=macoeffs:
        output['psi_'+str(s)]=vars()['psi_'+str(s)]
        s=s+1
    
    #might need to change this for no constant version
    
    
    #these G_s correspond to the matrices outlined at the top of page 337 of Hamilton
     
    
    s=1
    while s <= macoeffs:
        collector=np.zeros([numvars*numvars,numvars*(numvars*t+1)])
        i=1
        while i <= s:
            prodarray=np.zeros([numvars,1+numvars*t])
            p=1
            while p <=t:
                prodarray[:,1+(p-1)*numvars:1+p*numvars]=vars()['psi_'+str(s-i-p+1)].T
                p=p+1
                
            collector=collector+kron(vars()['psi_'+str(i-1)],prodarray)
            i=i+1
        vars()['G_'+str(s)]=collector
        s=s+1
    
    #output['collector']=collector
    
    #we should output the G_s to the output dictionary
    
    s=1
    while s<=macoeffs:
        output['G_'+str(s)]=vars()['G_'+str(s)]
        s=s+1
    
    
    #the Q_T matrix is necessary for calculation for either type of implied error bands
    
    Q_T=(1.0/(len(yt_)-len(betas)))*x_wcoll
    
    output['Q_T']=Q_T
    
    
    #now calculate the standard errors using equation from top of 337 of Hamilton
    #the (i,i) element of the standerrs_s matrix corresponds to the standard error 
    
    if orthog=='N' and errorbands=='Y':     
        s=1
        while s <=macoeffs:
            vars()['variances_'+str(s)]=(1.0/len(yt_))*dot(dot(vars()['G_'+str(s)],kron(omegahat,linalg.inv(Q_T))),vars()['G_'+str(s)].T)
            s=s+1
        
        #we should also output the standard errors
        
        s=1
        while s <=macoeffs:
            output['variances_'+str(s)]=vars()['variances_'+str(s)]
            s=s+1
        
        upci=np.zeros([numvars,nperiods])
        loci=np.zeros([numvars,nperiods])
        
        i=1
        while i <= macoeffs:
            upci[:,t+i]=1.96*sqrt(diag(vars()['variances_'+str(i)])[shockvar::numvars])
            loci[:,t+i]=-1.96*sqrt(diag(vars()['variances_'+str(i)])[shockvar::numvars])
            i=i+1
    
    
        output['upci']=upci
        output['loci']=loci
        
    
    ##Calculating standard errors requires a different method if the IRF has been orthogonalized
    
    ##vec() vs. vech() operator vech only includes elements on or below the main diagonal while vec includes all
    
    ##in order
    
    
    if orthog=='Y' and errorbands=='Y':
        
        upci=np.zeros([numvars,nperiods])
        loci=np.zeros([numvars,nperiods])
        
        
        #we will need the Dnplus matrix in order to calculate the standard errors with the orthogonalized VAR (Dnplus is the elimination matrix)
        #this requires the vech operator
        #this is the Duplication matrix
        #let us calculate the Elimination matrix first (going from vec to vech)
        
        #L_n=np.zeros([(1.0/2)*(numvars)*(numvars+1),numvars**2])#this is the elimination matrix
        
        
        #this method is taken from Magnus and Neudecker (1980)
        #########even though this technically is the elimination matrix, is NOT of the form specified in Hamilton
        #ident=np.eye(numvars)
        #j=0
        #while j < numvars:
        #    i=j
        #    while j <= i < numvars:
        #        u=np.zeros([(1.0/2)*numvars*(numvars+1),1])#initialize u_ij matrix
        #        u[(j)*numvars+(i+1)-(1.0/2)*(j+1)*j-1]=1 #because of 0 indexing, had to mess with formula from magnus (1980) paper
        #        L_n=L_n+kron(kron(u,ident[:,j].T),ident[:,i].T)
        #        i=i+1
        #    j=j+1
            
        #output['L_n']=L_n #this is the Elimination matrix
        
        #we can first determine the duplication matrix D_n
        D_nT=np.zeros([(1.0/2)*(numvars)*(numvars+1),numvars**2])

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
        
        output['D_n']=D_n #this is the duplication matrix pi=betas.ravel('F')[:,None]
        
        #this is the elimination matrix
        L_n=linalg.pinv(D_n)
        
        output['L_n']=L_n
        
        #in this case our B_0 matrix is simply equal to the inverse of our A matrix
        
        B_0=linalg.inv(A)
        
        B_0zeros=np.tril(B_0)
        
        output['B_0']=B_0
        
        output['B_0zeros']=B_0zeros
        
        #try this out
        #B_0=B_0zeros
        
        #we also need the S_B matrix, which solves the system vec(B_0.T)=dot(S_B.T,theta_B)+s_B
        #this method as currently coded is only meant for recursively identified shocks
        
        
        #i=0
        #c=0 #column counter
        #r=1
        #m=1
        #while i < numvars**2:
        #    r=m
        #    while r < numvars:
        #        S_B[r+i,c]=-1
        #        r=r+1 #this moves through the rows of the data
        #        c=c+1
        #    m=m+1 #this moves down the starting point for the first row to get a -1 value in its cth column every time we run through the numvars groups
        #    i=i+numvars
        
        
        #we can determine the S_B matrix and the S_BT matrix by manipulating the duplication matrix
        x=numvars #x keeps track of sum function to take out correct columns
        y=1 #y keeps track of lost columns
        collector=D_n
        collector=np.delete(collector,0,1)
        ###################here's the problem
        while x >= 2:
            collector=np.delete(collector,sum(b for b in range(x,numvars+1))-y,1)
            x=x-1
            y=y+1
        output['collector']=collector
        S_Bpre=collector
        
        output['S_Bpre']=S_Bpre
        
        
        #then for S_B we run through columns of collector matrix and set items following '1'ns in columns and set them to zero
        
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
        
        #turn all of the elements from 1 to -1
        S_BT=S_BT*(-1.0)
        S_B=S_B*(-1.0)
        
        output['S_B']=S_B
        output['S_BT']=S_BT
        
        #np.savetxt("G:\MacroProg\Output\S_B.txt",S_B,delimiter='\t')
        #np.savetxt("G:\MacroProg\Output\S_BT.txt",S_BT,delimiter='\t')
        
        
        S_D=np.zeros([numvars**2,numvars])
        
        i=0
        c=0 #column counter
        r=0
        while i < numvars**2:
            S_D[r+i,r]=1
            r=r+1 #this moves through the rows of the data
            i=i+numvars
            
        output['S_D']=S_D
        
        #need to intialize J matrix with dvech(omegahat)/dtheta_B and then append dvech(omegahat)/dtheta_D
        #remember that in this case A=linalg.inv(B_0)
        #########might want to change this for general case (not use A instead of lingalg.inv(B_0)
        J=-2.0*dot(dot(L_n,kron(omegahat,A)),S_B)
        J=append(J,dot(dot(L_n,kron(A,A)),S_D),axis=1)
        
        output['J']=J
        
        #in order to get G_B, G_B should be the first sum(x) rows of inv(J) matrix
        
        G_B=linalg.inv(J)[0:sum(x for x in range(0,numvars)),:]
        G_D=linalg.inv(J)[sum(x for x in range(0,numvars)):,:]
        
        output['G_B']=G_B
        output['G_D']=G_D
        
        
        #let us first compute the Ksi_pi and Ksi_sigma matrices for each s
        s=1
        while s <=macoeffs:
            vars()['Ksi_pi_'+str(s)]=dot(kron(np.eye(numvars),linalg.inv(B_0.T)),vars()['G_'+str(s)])
            vars()['Ksi_sigma_'+str(s)]=-dot(dot(kron(dot(vars()['psi_'+str(s)],A),linalg.inv(B_0.T)),S_BT),G_B)   #remember inv(B_0)=A
            s=s+1
            
        s=1
        while s<=macoeffs:
            output['Ksi_pi_'+str(s)]=vars()['Ksi_pi_'+str(s)]
            output['Ksi_sigma_'+str(s)]=vars()['Ksi_sigma_'+str(s)]
            s=s+1
           
        #now to compute the variance matrices
        
        s=1
        while s <=macoeffs:
            vars()['variances_'+str(s)]=(1.0/len(yt_))*(dot(dot(vars()['Ksi_pi_'+str(s)],kron(omegahat,linalg.inv(Q_T))),vars()['Ksi_pi_'+str(s)].T)+2.0*dot(dot(dot(dot(vars()['Ksi_sigma_'+str(s)],L_n),kron(omegahat,omegahat)),L_n.T),vars()['Ksi_sigma_'+str(s)].T))
            s=s+1
        
        
        i=1
        while i <= macoeffs:
            upci[:,t+i]=1.96*sqrt(diag(vars()['variances_'+str(i)])[shockvar::numvars])
            loci[:,t+i]=-1.96*sqrt(diag(vars()['variances_'+str(i)])[shockvar::numvars])
            i=i+1
    
    
        output['upci']=upci
        output['loci']=loci
        
        
        
        
        
    
    #######################now prepare the IRF that we choose for display
    
    #shockvar=int(raw_input("What variable would you like to add a shock to?: "))
    
    #shock[shockvar,t]=0.1
    
    solveirf=np.zeros([numvars,nperiods])
    
    #this follows the method outlined on pages 318 to 323
    
    if orthog=='Y':
        i=1
        while i<=macoeffs:
            solveirf[:,t+i]=dot(vars()['psi_'+str(i)],A[:,shockvar])
            i=i+1
        
    if orthog=='N':
        i=1
        while i<=macoeffs:
            solveirf[:,t+i]=vars()['psi_'+str(i)][:,shockvar]
            i=i+1
    
    
    output['solveirf']=solveirf
    
    if errorbands=='N':
        upci=np.zeros([numvars,nperiods])
        loci=np.zeros([numvars,nperiods])
    
    
    
    #showplot=raw_input("Would like you like to show the plot (Y/N)?: ")
    if showplot=="Y":
    	CEEplot.CEEplotter(solveirf,diff,loci,upci,names)
        
    return output
    

