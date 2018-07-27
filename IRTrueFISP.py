# Module to fit standard inversion recovery curve to True FISP

import numpy as np
import scipy.optimize

def SIeqn(paramsin, TIs, TR):
	# return signal intensity

	# paramsin is [T1,M0] - these are the unknowns
	# TR is TR in ms
	# TIs is array of TIs in ms

	T1=paramsin[0]
	M0=paramsin[1]

	# Calculate SI

	SI=M0*(1-2*np.exp(-1*TIs/T1)+np.exp((TIs-TR)/T1))
	return np.fabs(SI)

def fittingfun(TIs,TR,data):
	startguess=[500,10]
	bnds=((0,5000),(0,10000)) # Set upper and lower bounds for parameters
	fit=scipy.optimize.minimize(objfun,startguess,args=(TIs,TR,data),bounds=bnds, method='SLSQP',options={'ftol':1e-9,'disp':False,'eps':1e-10,'maxiter':1000})
	return fit


def objfun(paramsin,TIs,TR,data):
	# return chi**2
	SI=SIeqn(paramsin,TIs,TR)
	chi2=np.sum((data-SI)**2)
	return chi2