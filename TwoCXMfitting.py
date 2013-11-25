# Function to fit the TwoCXM model using scipy.optimize  
# Note that AIF is PLASMA

def TwoCXMfitting(t, AIF, uptake, toff):
	import numpy as np
	import scipy.optimize 
	import scipy.interpolate
	
	# Shift the AIF by the amount toff
	tnew = t - toff
	f = scipy.interpolate.PchipInterpolator(t,AIF)
	AIFnew = (t>toff)*f(t-toff)

	# Fit the TwoCXM model, stepping through vp
	vpmatrix=np.arange(0.01,1.01,0.01)
	
	# Parameters to fit are E, Fp, ve
	startguess=[0.5,0.5,0.5]  # Set starting guesses
	bnds=((0,1),(0,10),(0,2)) # Set upper and lower bounds for parameters
	resultsmatrix=np.zeros((len(vpmatrix),5))  # Initialise results array

	for i in range (0,len(vpmatrix)):
		Result=scipy.optimize.minimize(objfun,startguess,args=(np.array([vpmatrix[i]]),t,AIFnew,uptake),bounds=bnds, method='SLSQP',options={'ftol':1e-9,'disp':True,'eps':1e-10,'maxiter':1000})
		resultsmatrix[i,:]=(Result.x[0],Result.x[1],Result.x[2],vpmatrix[i],Result.fun)

	return resultsmatrix


def objfun(paramsin,vp,t,AIF,data):
    import numpy as np
    import TwoCXM
    allparams=np.concatenate((paramsin,vp))
    temp=data-TwoCXM.TwoCXM(allparams,t,AIF)
    return np.sqrt(sum(temp**2))