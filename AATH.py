# AATH module containing model function AATH and fitting functions for SI and Conc (AATHfitting...)
# Also includes tofts without vp to do toff estimation


def AATH(params,t,AIF, toff):
    import numpy as np
    import scipy.interpolate
    import matplotlib.pyplot as plt  
    
    # If toff value needs to be included (i.e. if not set to None), shift the AIF by the amount toff
    # Note when evaluating the model during fitting this shift has already been done
    if toff is not None:
        tnew = t - toff
        f=scipy.interpolate.interp1d(t,AIF,kind='cubic',bounds_error=False,fill_value=0)
        AIF = (t>=toff)*f(t-toff)

    #Test for trouble with fitting algorithm
    if np.isnan(np.sum(params)):
        F=np.zeros(len(AIF))
        return F
    # Assign the parameters to more meaningful names
    E, Fp, ve, vp = params
    Tc=vp/Fp
    
    #Calculate the IRF
    R=E*Fp*np.exp(-1*E*Fp*(t-Tc)/ve)
    if np.round(Tc/t[1])!=0:
        R[0:(np.round(Tc/t[1]))]=Fp
    #Calculate the convolution
    temp=np.convolve(AIF,R)*t[1]

    F=temp[0:len(t)]
    # plt.plot(AIF)
    # plt.plot(temp)
    # plt.plot(F)
    return F

# Model function for Kety with no vp to use for toff calculation
def Kety(params,t,AIF):
    import numpy as np
    import scipy.interpolate

    #test for trouble with fitting algorithm
    if np.isnan(np.sum(params)):
        F=np.zeros(len(AIF))
        return F
    # Assign parameter names
    #print(params)
    Ktrans, ve, toff = params

    # Shift the AIF by the amount toff
    tnew = t - toff
    f=scipy.interpolate.interp1d(t,AIF,kind='linear',bounds_error=False,fill_value=0)
    AIFnew = (t>toff)*f(t-toff)

    imp=Ktrans*np.exp(-1*Ktrans*t/ve); # Calculate the impulse response function
    convolution=np.convolve(AIFnew,imp) # Convolve impulse response with AIF
    G=convolution[0:len(t)]*t[1]
    return G


# Function to fit the AATH model using scipy.optimize  
# Note that AIF is PLASMA

# def AATHfittingSI(t, AIF, uptake, toff, baselinepts, TR, flip, T1base):
#     import numpy as np
#     import scipy.optimize 
#     import scipy.interpolate
#     import matplotlib.pyplot as plt 
#     import FLASH

#     #calculate M0 to use in conversion
#     rflip=flip*np.pi/180
#     R1base=1/T1base
#     base=np.mean(uptake[0:baselinepts])
#     print(base)
#     M0=base*(1-np.cos(rflip)*np.exp(-1*TR*R1base))/(np.sin(rflip)*(1-np.exp(-1*TR*R1base)))
#     print(M0)

#     # If toff is set to None, rather than a number, calculate it using Tofts without vp
#     plt.figure()

#     if toff is None:
#         firstthird=np.round(len(t)/3)
#         Ketystart=[0.5,0.5,t[1]] # set starting guesses for Ktrans, ve, toff
#         Ketybnds=((0,999),(0,2),(1E-05,1))
#         Ketyresult=scipy.optimize.minimize(Ketyobjfun,Ketystart,args=(t[0:firstthird],AIF[0:firstthird],uptake[0:firstthird],TR,flip,T1base,M0),bounds=Ketybnds,method='SLSQP')
#         toff=0
#         if not np.isnan(Ketyresult.x[2]):
#             toff=Ketyresult.x[2]
#         print(Ketyresult)
#         #concdata=Kety(Ketyresult.x,t,AIF)
#         #plt.plot(t,FLASH.Conc2SI(concdata,TR,flip,T1base,M0),'b')

#     if toff==0:
#         AIFnew=AIF

#     else:
#         # Shift the AIF by the amount toff
#         tnew = t - toff
#         f=scipy.interpolate.interp1d(t,AIF,kind='cubic',bounds_error=False,fill_value=0)
#         AIFnew = (t>=toff)*f(t-toff)

#     # Fit the AATH model, stepping through vp
#     vpmatrix=np.arange(0.01,1,0.01)
#     #vpmatrix=np.array([0.4])    
    
#     # Parameters to fit are E, Fp, ve
#     startguess=[0.1,0.1,0.1]  # Set starting guesses
#     bnds=((1e-5,1),(1e-05,10),(1e-05,1)) # Set upper and lower bounds for parameters
#     resultsmatrix=np.zeros((len(vpmatrix),6))  # Initialise results array

#     for i in range (0,len(vpmatrix)):
#         Result=scipy.optimize.minimize(objfun,startguess,args=(np.array([vpmatrix[i]]),t,AIFnew,uptake,TR,flip,T1base,M0),bounds=bnds, method='SLSQP',options={'ftol':1e-14,'disp':False, 'eps':1e-14, 'maxiter':1000})
#         resultsmatrix[i,:]=(Result.x[0],Result.x[1],Result.x[2],vpmatrix[i],Result.fun,toff)

#     #print(resultsmatrix)
#     bestindex=np.nanargmin(resultsmatrix[:,4])
#     bestresult=resultsmatrix[bestindex,:]

#     #plt.plot(resultsmatrix[:,4])
#     #plt.figure()
#     plt.plot(t,uptake,'x')
#     plt.plot(t,FLASH.Conc2SI(AATH(bestresult[0:4],t,AIF,toff),TR,flip,T1base,M0),'r')

#     return bestresult

def AATHfittingConc(t, AIF, uptake, toff):
    import numpy as np
    import scipy.optimize 
    import scipy.interpolate
    import matplotlib.pyplot as plt

    # If toff is set to None, rather than a number, calculate it using Tofts without vp from the first third of the curve
    #plt.figure()

    firstthird=np.round(len(t)/3)
    if toff is None:
        Ketystart=[0.05,0.5,t[1]] # set starting guesses for Ktrans, ve, toff
        Ketybnds=((0.00001,999),(0.00001,2),(0.00001,20))
        Ketyresult=scipy.optimize.minimize(KetyobjfunConc,Ketystart,args=(t[0:firstthird],AIF[0:firstthird],uptake[0:firstthird]),bounds=Ketybnds,method='SLSQP')
        toff=0
        if not np.isnan(Ketyresult.x[2]):
            toff=Ketyresult.x[2]
        #plt.plot(t,Kety(Ketyresult.x[0:4],t,AIF))
        #print(Ketyresult.x)
        print('Success? '+str(Ketyresult.success)+' toff='+str(Ketyresult.x[2]))
    
    # Shift the AIF by the amount toff
    tnew = t - toff
    f=scipy.interpolate.interp1d(t,AIF,kind='linear',bounds_error=False,fill_value=0)
    AIFnew = (t>=toff)*f(t-toff)

    # Fit the AATH model, stepping through vp
    vpmatrix=np.arange(0.01,1,0.01) #was 0.01 start
    # Parameters to fit are E, Fp, ve
    startguess=[0.5,0.5,0.5]  # Set starting guesses
    bnds=((0.00001,1),(0.00001,10),(0.00001,3)) # Set upper and lower bounds for parameters
    resultsmatrix=np.zeros((len(vpmatrix),6))  # Initialise results array

    for i in range (0,len(vpmatrix)):
        Result=scipy.optimize.minimize(objfunConc,startguess,args=(np.array([vpmatrix[i]]),t,AIFnew,uptake),bounds=bnds, method='SLSQP',options={'ftol':1e-14,'disp':False,'eps':1e-14,'maxiter':1000})
        resultsmatrix[i,:]=(Result.x[0],Result.x[1],Result.x[2],vpmatrix[i],Result.fun,toff)
    
    #print(resultsmatrix)
    bestindex=np.nanargmin(resultsmatrix[:,4])
    bestresult=resultsmatrix[bestindex,:]
    print(bestresult)
    #plt.plot(resultsmatrix[:,4])
    #plt.figure()
    #plt.plot(t,uptake,'x')
    #plt.plot(t,AATH(bestresult[0:4],t,AIF,toff))
    
    return bestresult

def Ketyobjfun(paramsin,t,AIF,data,TR,flip,T1base,M0):
    import numpy as np
    import FLASH
    concdata=Kety(paramsin,t,AIF)
    temp=data-FLASH.Conc2SI(concdata,TR,flip,T1base,M0)
    return np.sqrt(np.sum(temp**2))

def KetyobjfunConc(paramsin,t,AIF,data):
    import numpy as np
    temp=data-Kety(paramsin,t,AIF)
    return np.sqrt(np.sum(temp**2))

def objfun(paramsin,vp,t,AIF,data,TR,flip,T1base,M0):
    import numpy as np
    import FLASH
    allparams=np.concatenate((paramsin,vp))
    concdata=AATH(allparams,t,AIF,None)
    temp=data-FLASH.Conc2SI(concdata,TR,flip,T1base,M0)
    return np.sqrt(np.sum(temp**2))
    
def objfunConc(paramsin,vp,t,AIF,data):
    import numpy as np
    allparams=np.concatenate((paramsin,vp))
    temp=data-AATH(allparams,t,AIF,None)
    return np.sqrt(np.sum(temp**2))
