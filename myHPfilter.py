def myHPfilter(y,w):
    #Based on Ivailo Izvorski's Matlab Code
    #see http://dge.repec.org/codes.html

    # w - smoothing parameter
    # y - original series
    # s - filtered series

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
