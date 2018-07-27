# FLASH signal equation, conversion from SI to concentration and vice-versa

import numpy as np
import scipy.optimize

def SIeqn(paramsin, flip, TR):
	M0=paramsin[0]
	T1=paramsin[1]
	rflip=flip*np.pi/180
	SI=M0*np.sin(rflip)*(1-np.exp(-TR/T1))/(1-np.cos(rflip)*np.exp(-TR/T1))
	return SI


def fittingfun(flips,TR,data,startguess=[10000,500]):
	#startguess is [M0,T1]
	bnds=((0,1E10),(0,100000)) # Set upper and lower bounds for parameters
	fit=scipy.optimize.minimize(objfun,startguess,args=(flips,TR,data),bounds=bnds, method='SLSQP',options={'ftol':1e-9,'disp':False,'maxiter':100000})
	#fit=scipy.optimize.minimize(objfun,startguess,args=(flips,TR,data),bounds=bnds, method='Nelder-Mead',options={'ftol':1e-9,'disp':False,'eps':1e-10,'maxiter':100000})
	return fit


def objfun(paramsin,flips,TR,data):
	# return chi**2
	SI=SIeqn(paramsin,flips,TR)
	chi2=np.sum((data-SI)**2)
	return chi2


# Conversion from SI to concentration
def SI2Conc(SIcurve,TR,flip,T1base,baselinepts,M0,use2ndT1=0):
	# TR and T1 should be in s, resulting delta R1 will be in s^-1
	import numpy as np

	# Convert flip angle to radians
	rflip=flip*np.pi/180
	# Convert T1 to R1
	R1base=1/T1base

	# If M0 isn't specified, calculate from baseline (missing point 1)
	if M0 is None:
		base=np.mean(SIcurve[1:baselinepts])
		M0=base*(1-np.cos(rflip)*np.exp(-1*TR*R1base))/(np.sin(rflip)*(1-np.exp(-1*TR*R1base)))
	
	# Now calculate the R1 curve
	R1=np.log(((M0*np.sin(rflip))-SIcurve)/(M0*np.sin(rflip)-(SIcurve*np.cos(rflip))))*(-1/TR)

	# And finally the delta R1 curve
	if use2ndT1==1:
		R1base=np.mean(R1[1:baselinepts]) #If using the post-contrast T1 measurement, subtract a suitable R1 to make the baseline=0
	H=R1-R1base
	#return H
	return H

def Conc2SI(deltaR1,TR,flip,T1base,M0,use2ndT1=0):
	# TR and T1 should be in s, delta R1 in s^-1
	import numpy as np

	# Convert flip angle to radians
	rflip=flip*np.pi/180
	# Convert T1 base to R1 base
	R1base=1/T1base

	# Convert deltaR1 curve to R1 curve
	if use2ndT1==1:
		R1base=R1base-np.mean(deltaR1[-15:]) #If using second T1, calculate the necessary R1 to be added to delta R1 to match R1 at the end (confusingly R1base in this case)
	R1curve=deltaR1+R1base
	# Convert to SI
	SI=M0*np.sin(rflip)*(1-np.exp(-TR*R1curve))/(1-np.cos(rflip)*np.exp(-TR*R1curve))
	return SI

def CalcM0(SI, TR, flip, T1):
	# TR and T1 should be in s
	import numpy as np
	#convert flip angle to radians
	rflip=flip*np.pi/180
	#Calculate M0
	M0=SI*(1-np.cos(rflip)*np.exp(-TR/T1))/(np.sin(rflip)*(1-np.exp(-TR/T1)))
	return M0