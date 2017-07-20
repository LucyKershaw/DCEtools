# Model function for fitting extended Kety in Concentration space

# Model function for extended Kety (with vp)
def ExtKety(params,t,AIF):
    import numpy as np
    import scipy.interpolate
    #Test for trouble with fitting algorithm
    if np.isnan(np.sum(params)):
        F=np.zeros(len(AIF))
        return F

    # Assign parameter names
    Ktrans, ve, vp, toff = params

    # Shift the AIF by the amount toff
    tnew = t - toff
    f=scipy.interpolate.interp1d(t,AIF,kind='linear',bounds_error=False,fill_value=0)
    AIFnew = (t>toff)*f(t-toff)

    imp=Ktrans*np.exp(-1*Ktrans*t/ve); # Calculate the impulse response function
    convolution=np.convolve(AIFnew,imp) # Convolve impulse response with AIF
    G=convolution[0:len(t)]*t[1]+vp*AIFnew  #adjust for temporal resolution and add Cpvp term

    return G

def ExtKetyobjfun(paramsin,t,AIF,data):
    import numpy as np
    temp=data-ExtKety(paramsin,t,AIF)
    return np.sqrt(np.sum(temp**2))

def ExtKetyobjfunSI(paramsin,t,AIF,data, TR, flip, T1base, M0):
    import numpy as np
    import FLASH
    concdata=ExtKety(paramsin,t,AIF)
    temp=data-FLASH.Conc2SI(concdata,TR,flip,T1base,M0)
    return np.sqrt(np.sum(temp**2))


def ExtKetyfittingConc(t, AIF, uptake):
    import numpy as np
    import scipy.optimize 
    import scipy.interpolate
    import matplotlib.pyplot as plt

    #plt.figure()

    Ketystart=np.array((0.2,0.1,0.01,t[1])) # set starting guesses for Ktrans, ve, vp, toff
    Ketybnds=((0.00001,999),(0.00001,1),(0.00001,1),(0.00001,20))
    ExtKetyresult=scipy.optimize.minimize(ExtKetyobjfun,Ketystart,args=(t,AIF,uptake),bounds=Ketybnds,method='SLSQP',options={'disp':False})
            
    #plt.plot(t,ExtKety(ExtKetyresult.x[0:4],t,AIF))
    #plt.plot(t,uptake)
    #print(ExtKetyresult.x)
    print(ExtKetyresult.success)
    return(ExtKetyresult.x)

def ExtKetyfittingSI(t, AIF, uptake, baselinepts, TR, flip, T1base):
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

    Ketystart=np.array((0.1,0.01,0.01,t[1])) # set starting guesses for Ktrans, ve, vp, toff
    Ketybnds=((0.00001,999),(0.00001,1),(0.00001,1),(0.00001,30))
    ExtKetyresult=scipy.optimize.minimize(ExtKetyobjfunSI,Ketystart,args=(t,AIF,uptake,TR,flip,T1base,M0),bounds=Ketybnds,method='SLSQP',options={'disp':False})

    print(ExtKetyresult.success)
    return(ExtKetyresult.x)



