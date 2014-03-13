
# Module to deal with fitting the TH model, in the fourier domain
import numpy as np
import matplotlib.pyplot as plt
#Model function
def TH(params,t,AIF,toff):
	#parameters are [E,Fp,ve,vp] as usual
	E=params[0]
	Fp=params[1]
	ve=params[2]
	vp=params[3]

	#First calculate model parameters alpha, beta, chi

	#alpha=PS/Fp - plug flow therefore Renkin-Crone (E=1-e^(-PS/Fp)) applicable.
	alpha=-1*np.log(1-E)
	beta=vp/Fp
	chi=ve/vp

	#calculate zero padding 
	N=len(t)
	tres=t[1]
	M=int(N+(6*ve/(E*Fp*tres)))

	#construct frequency vector from deltaf=1/Mdeltat
	deltaf=1/(M*tres)
	f=deltaf*np.arange(0,M,1)

	#Now calculate H, the impulse response in the Fourier domain H=H(2*pi*i*deltaf*m)
	s=2*np.pi*1j*f
	B1=(beta*s)+alpha
	H=((np.exp(-1*B1)-1)*B1*((alpha*ve)+vp*((s*chi*beta)+alpha)))/(alpha*alpha*(np.exp(-1*B1)-1)-(B1*s*beta*(chi*B1)+alpha))
	plt.plot(H)

	#need the Fourier transform of the AIF
	AIFT=np.fft.fft(AIF)

	#Multiply by H
	CtT=AIFT*H
	plt.plot(CtT)

	#Transform back
	Ct=np.fft.ifft(CtT,M)

	return Ct

#fitting function in conc space
def THfittingConc(t,AIF,uptake,toff):
	pass
#Objective function
def THobjfunConc(params,t,AIF,data):
	pass
