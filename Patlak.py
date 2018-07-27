# Function for calculating and fitting Patlak model

import matplotlib.pyplot as plt
import scipy.integrate

# Model function for Patlak

def Patlak(params,t,AIF,plot=0):
    import numpy as np
    import scipy.interpolate
    #Test for trouble with fitting algorithm
    if np.isnan(np.sum(params)):
        F=np.zeros(len(AIF))
        return F

    # Assign parameter names
    Ktrans, vp = params

    # Shift the AIF by the amount toff
    toff=0
    tnew = t - toff
    f=scipy.interpolate.interp1d(t,AIF,kind='linear',bounds_error=False,fill_value=0)
    AIFnew = (t>toff)*f(t-toff)

    #convolution=np.convolve(AIFnew,Ktrans) # Convolve impulse response with AIF
    convolution=Ktrans*scipy.integrate.cumtrapz(AIFnew,x=tnew,initial=0)
    #print(convolution)
    #G=convolution[0:len(t)]*t[1]+(vp*AIFnew)  #adjust for temporal resolution and add Cpvp term
    G=convolution+(AIFnew*vp)

    if plot==1:
        plt.plot(t,G,'k')
        plt.plot(t,convolution,'b')
        plt.plot(t,vp*AIFnew,'r')

    return G

def Patlakobjfun(paramsin,t,AIF,data):
    import numpy as np
    temp=data-Patlak(paramsin,t,AIF)
    return np.sqrt(np.sum(temp**2))

def PatlakobjfunSI(paramsin,t,AIF,data, TR, flip, T1base, M0):
    import numpy as np
    import FLASH
    concdata=Patlak(paramsin,t,AIF)
    temp=data-FLASH.Conc2SI(concdata,TR,flip,T1base,M0)
    return np.sqrt(np.sum(temp**2))


def PatlakfittingConc(t, AIF, uptake):
    import numpy as np
    import scipy.optimize 
    import scipy.interpolate
    import matplotlib.pyplot as plt

    #plt.figure()

    Patlakstart=np.array((0.3,0.4)) # set starting guesses for Ktrans, vp, toff
    Patlakbnds=((1E-10,999),(1E-10,1))
    Patlakresult=scipy.optimize.minimize(Patlakobjfun,Patlakstart,args=(t,AIF,uptake),bounds=Patlakbnds,method='SLSQP',options={'ftol':1e-12,'disp':3,'maxiter':100000})
            
    #plt.plot(t,Patlak(Patlakresult.x[0:4],t,AIF))
    #plt.plot(t,uptake)
    #print(Patlaresult.x)
    print(Patlakresult)
    return(Patlakresult.x)

def PatlakfittingSI(t, AIF, uptake, baselinepts, TR, flip, T1base):
    import numpy as np
    import scipy.optimize 
    import scipy.interpolate
    import matplotlib.pyplot as plt 
    import FLASH

    #calculate M0 to use in conversion
    rflip=flip*np.pi/180
    R1base=1/T1base
    base=np.mean(uptake[0:baselinepts])
    #print(base)
    M0=base*(1-np.cos(rflip)*np.exp(-1*TR*R1base))/(np.sin(rflip)*(1-np.exp(-1*TR*R1base)))
    #print(M0)

    Patlakstart=np.array((0.1,0.01,t[1])) # set starting guesses for Ktrans, vp, toff
    Patlakbnds=((0.00001,999),(0.00001,1),(0.00001,30))
    Patlakresult=scipy.optimize.minimize(PatlakobjfunSI,Patlakstart,args=(t,AIF,uptake,TR,flip,T1base,M0),bounds=Patlakbnds,method='SLSQP',options={'disp':False})

    print(Patlakresult.success)
    return(Patlakresult.x)