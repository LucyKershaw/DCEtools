# Module to define functions necessary to model an SRTurboFLASH sequence

import numpy as np
import scipy.optimize

def SIeqn(paramsin,TIs,TR,flip,n):
	#return signal intensity

	#paramsin is [T1,M0] - these are the unknowns
	#n is steps to central k space
	#TR is the echospacing TR in ms
	#TIs is array of TIs in ms
	
	TIscorr=TIs-(n*TR) #TI in signal equation is to first alpha pulse, Siemens reports to central k space
	radflip=flip*np.pi/180
	T1=paramsin[0]
	M0=paramsin[1]
	ETR=np.exp(-1*TR/T1)
	ETI=np.exp(-1*TIscorr/T1)
	cosETR=ETR*np.cos(radflip)
	SI=M0*np.sin(radflip)*((1-ETI)*(cosETR**(n-1))+(1-ETR)*(1-(cosETR**(n-1)))/(1-cosETR))
	return SI

def fittingfun(TIs,TR,flip,n,data):
	startguess=[1000,1500]
	bnds=((0,3000),(0,1000000)) # Set upper and lower bounds for parameters
	fit=scipy.optimize.minimize(objfun,startguess,args=(TIs,TR,flip,n,data),bounds=bnds, method='SLSQP',options={'ftol':1e-14,'disp':True,'eps':1e-10,'maxiter':1000})
	return fit


def objfun(paramsin,TIs,TR,flip,n,data):
	# return chi**2
	SI=SIeqn(paramsin,TIs,TR,flip,n)
	chi2=np.sum((data-SI)**2)
	return chi2

