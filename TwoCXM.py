# TwoCXM module containing model function TwoCXM and fitting function TwoCXMfitting
# Also includes tofts without vp to do toff estimation


def TwoCXM(params,t,AIF, toff):
    import numpy as np
    import scipy.interpolate
    from scipy.signal import fftconvolve as conv
    import matplotlib.pyplot as plt    
    
    # If toff value needs to be included (i.e. if not set to None), shift the AIF by the amount toff
    if toff is not None:
        tnew = t - toff
        f = scipy.interpolate.PchipInterpolator(t,AIF)
        AIF = (t>=toff)*f(t-toff)

    # Assign the parameters to more meaningful names
    E, Fp, ve, vp = params

    #First calculate the parameters TB, TE and Tp
    Tp=(vp/Fp)*(1-E)
    TE=ve*(1-E)/(E*Fp)
    TB=vp/Fp

    #And then the IRF parameters A, Kplus, Kminus
    Kplus=0.5*((1/Tp) + (1/TE) + np.sqrt(((1/Tp) + (1/TE))**2 - (4/(TE*TB))))
    Kminus=0.5*((1/Tp) + (1/TE) - np.sqrt(((1/Tp) + (1/TE))**2 - (4/(TE*TB))))
    A=(Kplus - (1/TB))/(Kplus - Kminus)

    #Calculate the IRF
    expKplus=np.exp(-t*Kplus)
    expKminus=np.exp(-t*Kminus)

    R=expKplus + A*(expKminus - expKplus)
    #Calculate the convolution
    temp=conv(AIF,R)*t[1]
    F=Fp*temp[0:len(t)]
    # plt.plot(AIF)
    # plt.plot(temp)
    # plt.plot(F)
    return F

# Model function for Kety with no vp to use for toff calculation
def Kety(params,t,AIF):
    import numpy as np
    import scipy.interpolate
    from scipy.signal import fftconvolve as conv

    # Assign parameter names
    Ktrans, ve, toff = params

    # Shift the AIF by the amount toff
    tnew = t - toff
    f = scipy.interpolate.PchipInterpolator(t,AIF)
    AIFnew = (t>toff)*f(t-toff)

    imp=ve*np.exp(-1*Ktrans*t/ve); # Calculate the impulse response function
    convolution=conv(AIFnew,imp) # Convolve impulse response with AIF
    G=convolution[0:len(t)]*t[1]
    return G


# Function to fit the TwoCXM model using scipy.optimize  
# Note that AIF is PLASMA

def TwoCXMfittingSI(t, AIF, uptake, toff, baselinepts, TR, flip, T1base):
    import numpy as np
    import scipy.optimize 
    import scipy.interpolate

    #calculate M0 to use in conversion
    rflip=flip*np.pi/180
    R1base=1/T1base
    base=np.mean(uptake[0:baselinepts])
    print(base)
    M0=base*(1-np.cos(rflip)*np.exp(-1*TR*R1base))/(np.sin(rflip)*(1-np.exp(-1*TR*R1base)))
    print(M0)
    print('This code')

    # If toff is set to None, rather than a number, calculate it using Tofts without vp

    if toff is None:
        Ketystart=[0.5,0.5,t[1]] # set starting guesses for Ktrans, ve, toff
        Ketybnds=((0,999),(0,2),(0,1))
        Ketyresult=scipy.optimize.minimize(Ketyobjfun,Ketystart,args=(t,AIF,uptake,TR,flip,T1base,M0),bounds=Ketybnds,method='SLSQP')
        toff=Ketyresult.x[2]
        print(toff)
    
    if toff==0:
        AIFnew=AIF

    else:
        # Shift the AIF by the amount toff
        tnew = t - toff
        f = scipy.interpolate.PchipInterpolator(t,AIF)
        AIFnew = (t>=toff)*f(t-toff)

    # Fit the TwoCXM model, stepping through vp
    vpmatrix=np.arange(0.01,0.2,0.001)
    #vpmatrix=np.array([0.4])    
    
    # Parameters to fit are E, Fp, ve
    startguess=[0.1,0.1,0.1]  # Set starting guesses
    bnds=((1e-5,1),(1e-05,10),(1e-05,1)) # Set upper and lower bounds for parameters
    resultsmatrix=np.zeros((len(vpmatrix),5))  # Initialise results array

    for i in range (0,len(vpmatrix)):
        Result=scipy.optimize.minimize(objfun,startguess,args=(np.array([vpmatrix[i]]),t,AIFnew,uptake,TR,flip,T1base,M0),bounds=bnds, method='SLSQP',options={'ftol':1e-14,'disp':False,'maxiter':1000})
        resultsmatrix[i,:]=(Result.x[0],Result.x[1],Result.x[2],vpmatrix[i],Result.fun)

    print(resultsmatrix)
    bestindex=np.nanargmin(resultsmatrix[:,4])
    bestresult=resultsmatrix[bestindex,:]
    return bestresult

def TwoCXMfittingConc(t, AIF, uptake, toff):
    import numpy as np
    import scipy.optimize 
    import scipy.interpolate
    import matplotlib.pyplot as plt

    # If toff is set to None, rather than a number, calculate it using Tofts without vp form the first third of the curve

    firstthird=np.round(len(t)/3)
    if toff is None:
        Ketystart=[0.5,0.5,t[1]] # set starting guesses for Ktrans, ve, toff
        Ketybnds=((0,999),(0,2),(1E-5,1))
        Ketyresult=scipy.optimize.minimize(KetyobjfunConc,Ketystart,args=(t[0:firstthird],AIF[0:firstthird],uptake[0:firstthird]),bounds=Ketybnds,method='SLSQP')
        toff=Ketyresult.x[2]
        print(Ketyresult)
    
    # Shift the AIF by the amount toff
    tnew = t - toff
    f = scipy.interpolate.PchipInterpolator(t,AIF)
    AIFnew = (t>=toff)*f(t-toff)

    # Fit the TwoCXM model, stepping through vp
    vpmatrix=np.arange(0.3,0.5,0.01) #was 0.01 start
    # Parameters to fit are E, Fp, ve
    startguess=[0.5,0.5,0.5]  # Set starting guesses
    bnds=((1e-5,1),(1e-5,10),(1e-5,3)) # Set upper and lower bounds for parameters
    resultsmatrix=np.zeros((len(vpmatrix),5))  # Initialise results array

    for i in range (0,len(vpmatrix)):
        Result=scipy.optimize.minimize(objfunConc,startguess,args=(np.array([vpmatrix[i]]),t,AIFnew,uptake),bounds=bnds, method='SLSQP',options={'ftol':1e-14,'disp':True,'eps':1e-14,'maxiter':1000})
        resultsmatrix[i,:]=(Result.x[0],Result.x[1],Result.x[2],vpmatrix[i],Result.fun)
    
    print(resultsmatrix)
    bestindex=np.nanargmin(resultsmatrix[:,4])
    finalresult=np.concatenate((resultsmatrix[bestindex,:],[toff]))
    print(finalresult)
    plt.plot(resultsmatrix[:,4])
    plt.figure()
    plt.plot(uptake,'x')
    plt.plot(TwoCXM(finalresult[0:4],t,AIF,toff))
    if toff is None:
        plt.plot(Kety(Ketyresult.x[0:4],t,AIF))

    return finalresult

def Ketyobjfun(paramsin,t,AIF,data,TR,flip,T1base,M0):
    import numpy as np
    import FLASH
    concdata=Kety(paramsin,t,AIF)
    temp=data-FLASH.Conc2SI(concdata,TR,flip,T1base,M0)
    return np.sqrt(sum(temp**2))

def KetyobjfunConc(paramsin,t,AIF,data):
    import numpy as np
    temp=data-Kety(paramsin,t,AIF)
    return np.sqrt(sum(temp**2))

def objfun(paramsin,vp,t,AIF,data,TR,flip,T1base,M0):
    import numpy as np
    import FLASH
    allparams=np.concatenate((paramsin,vp))
    concdata=TwoCXM(allparams,t,AIF,None)
    temp=data-FLASH.Conc2SI(concdata,TR,flip,T1base,M0)
    return np.sqrt(sum(temp**2))
    
def objfunConc(paramsin,vp,t,AIF,data):
    import numpy as np
    allparams=np.concatenate((paramsin,vp))
    #temp=data-TwoCXM(allparams,t,AIF,None)
    temp=TwoCXM(allparams,t,AIF,None)
    #return np.median(np.absolute(temp))
    #return np.sqrt(sum(temp**2))
    return -1*sum((data[10:]*np.log(temp[10:]))-data[10:])
