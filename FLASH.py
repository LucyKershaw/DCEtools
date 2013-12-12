# FLASH signal equation - conversion from SI to concentration and vice-versa

def SIeqn(M0, flip, TR, T1base, ):
	pass




# Conversion from SI to concentration
def SI2Conc(SIcurve,TR,flip,T1base,baselinepts,M0):
	# TR and T1 should be in s, resulting delta R1 will be in s^-1
	import numpy as np

	# Convert flip angle to radians
	rflip=flip*np.pi/180
	# Convert T1 to R1
	R1base=1/T1base

	# If M0 isn't specified, calculate from baseline
	if M0 is None:
		base=np.mean(SIcurve[0:baselinepts])
		M0=base*(1-np.cos(rflip)*np.exp(-1*TR*R1base))/(np.sin(rflip)*(1-np.exp(-1*TR*R1base)))
	
	# Now calculate the R1 curve
	R1=np.log(((M0*np.sin(rflip))-SIcurve)/(M0*np.sin(rflip)-(SIcurve*np.cos(rflip))))*(-1/TR)

	# And finally the delta R1 curve
	H=R1-R1base
	return H


def Conc2SI(deltaR1,TR,flip,T1base,M0):
	# TR and T1 should be in s, delta R1 in s^-1
	import numpy as np

	# Convert flip angle to radians
	rflip=flip*np.pi/180
	# Convert T1 base to R1 base and TR to seconds
	R1base=1/T1base

	# Convert deltaR1 curve to R1 curve
	R1curve=deltaR1+R1base
	# Convert to SI
	SI=M0*np.sin(rflip)*(1-np.exp(-TR*R1curve))/(1-np.cos(rflip)*np.exp(-TR*R1curve))
	return SI