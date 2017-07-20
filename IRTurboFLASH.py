# Module to define functions necessary to model an IRTurboFLASH sequence

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

def SIeqn(paramsin,TIs,sequenceparams,abs=1):
	# return signal intensity

	# paramsin is [T1,M0] - these are the unknowns
	# Sequenceparams are [flipangle, n, TR, N, bigTR]
	# n is steps to central k space
	# TR is TR in ms
	# N is total k-space steps
	# bigTR is the shot interval
	# TIs is array of TIs in ms

	flip, n, TR, N, bigTR = sequenceparams 
	
	TIscorr=TIs-(n*TR) #TI in signal equation is to first alpha pulse, Siemens and Philips report to central k space
	radflip=flip*np.pi/180
	TD=bigTR-TIscorr-(N*TR)

	T1=paramsin[0]
	M0=paramsin[1]
	
	# Calculate constants to use later
	a=(np.exp(-1*(TR/T1)))*np.cos(radflip)
	b=1-np.exp(-1*TR/T1)
	C=(a**(N-1))*(1-np.exp(-1*TIscorr/T1))
	A=(1-a**(N-1))/(1-a)
	D=(np.cos(radflip))*np.exp(-1*TD/T1)
	B=(a**(N-1))*np.exp(-1*TIscorr/T1)

	E=C+(b*A)-(1/np.cos(radflip))

	firstfraction=-1*(((D*E)+1)/(1+(B*D)))

	#SI=np.fabs((M0*np.sin(radflip)*((firstfraction*(np.exp(-1*TIscorr/T1)*(a**(n-1))+((1-np.exp(-1*TIscorr/T1)*(a**(n-1))+(b*(1-(a**(n-1))/(1-a)))))))
	if abs==1:
		SI=np.fabs((M0*np.sin(radflip))*((firstfraction*np.exp(-1*TIscorr/T1)*a**(n-1))+((1-np.exp(-1*TIscorr/T1))*a**(n-1))+(b*(1-a**(n-1))/(1-a))))
	else:
		SI=((M0*np.sin(radflip))*((firstfraction*np.exp(-1*TIscorr/T1)*a**(n-1))+((1-np.exp(-1*TIscorr/T1))*a**(n-1))+(b*(1-a**(n-1))/(1-a))))
	return SI

def fittingfun(TIs,sequenceparams,data):
	startguess=[1000,10]
	bnds=((0,3000),(0,10000)) # Set upper and lower bounds for parameters
	fit=scipy.optimize.minimize(objfun,startguess,args=(TIs,sequenceparams,data),bounds=bnds, method='SLSQP',options={'ftol':1e-9,'disp':False,'eps':1e-10,'maxiter':1000})
	return fit


def objfun(paramsin,TIs,sequenceparams,data):
	# return chi**2
	SI=SIeqn(paramsin,TIs,sequenceparams)
	chi2=np.sum((data-SI)**2)
	return chi2


