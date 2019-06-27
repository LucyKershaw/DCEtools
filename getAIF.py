# Module to do AIF choosing from pixels

def getAIF(dynimages, slice, TR, flip, peakframe, numvessels=2, blinepts=15):
	# TR should be in seconds
	import numpy as np
	import matplotlib.pyplot as plt
	#import pyqtgraph as pg
	import FLASH
	
	# Show image and click the vessels
	h1=plt.figure()
	plt.imshow(dynimages[:,:,slice,peakframe])
	vesselpos=plt.ginput(n=numvessels,timeout=100, show_clicks=True,mouse_add=1, mouse_pop=3, mouse_stop=2) # these are not integers
	plt.close(h1)
	AIF=np.zeros((numvessels,dynimages.shape[3]))

	# For each vessel location....
	for v in range(0,numvessels):
		# Display a 9x9 grid of AIFs from around the points
		plt.ioff()
		g,axarr=plt.subplots(9,9,sharex=False,sharey=True,figsize=[15,10],tight_layout='True')
		

		#j is columns on the plot, i is rows

		for i in range(-4,5):
			for j in range(-4,5):
				vesselx=np.int(np.round(vesselpos[v][0]))
				vessely=np.int(np.round(vesselpos[v][1]))
				axarr[i+4,j+4].plot(dynimages[vessely+i,vesselx+j,slice,:])
				axarr[i+4,j+4].set_xticklabels([])
				axarr[i+4,j+4].tick_params(labelsize=8)
				axarr[i+4,j+4].set_xbound([0,dynimages.shape[3]])
				axarr[i+4,j+4].text(0.9,0.7,str(9*(i+4)+j+4),horizontalalignment='right',transform=axarr[i+4,j+4].transAxes)
		#g.draw()
		g.show()

		# Convert to concentration and plot again
		h,axarr=plt.subplots(9,9,sharex=False,sharey=True,figsize=[15,10],tight_layout='True')
		
		for i in range(-4,5):
			for j in range(-4,5):
				vesselx=np.int(np.round(vesselpos[v][0]))
				vessely=np.int(np.round(vesselpos[v][1]))
				axarr[i+4,j+4].plot(FLASH.SI2Conc(dynimages[vessely+i,vesselx+j,slice,:],TR,flip,1.68,blinepts,None))
				axarr[i+4,j+4].set_xticklabels([])
				axarr[i+4,j+4].tick_params(labelsize=8)
				axarr[i+4,j+4].set_xbound([0,dynimages.shape[3]])
				#axarr[i+4,j+4].set_ybound([-1,20])
				axarr[i+4,j+4].text(0.9,0.7,str(9*(i+4)+j+4),horizontalalignment='right',transform=axarr[i+4,j+4].transAxes)
		ysize=axarr[4,4].get_ybound()
		for i in range(-4,5):
			for j in range(-4,5):
				axarr[i+4,j+4].set_ybound(ysize)

		#g.draw()
		h.show()
		plt.ion()

		# Now choose which subplots to use
		ids=input("Input the id number for the chosen AIF plots, separated by commas: ")
		ids=ids.split(',')
		ChosenAIFs=[int(item) for item in ids]
		plt.close(h)
		plt.close(g)

		# Work out i,js
		AIFij=[[int((x-np.remainder(x,9))/9),np.remainder(x,9)] for x in ChosenAIFs]
		#print(AIFij)
		chosenSIcurves=np.zeros((len(ChosenAIFs),dynimages.shape[3]))
		chosenConccurves=np.zeros((len(ChosenAIFs),dynimages.shape[3]))

		# And extract them
		for p in range(0,len(AIFij)):
			chosenSIcurves[p,:]=np.squeeze(dynimages[vessely+AIFij[p][0]-4,vesselx+AIFij[p][1]-4,slice,:])
			chosenConccurves[p,:]=FLASH.SI2Conc(chosenSIcurves[p,:],TR,flip,1.68,blinepts,None)
			print(str(slice)+', '+str(vessely+AIFij[p][0]-4)+', '+str(vesselx+AIFij[p][1]-4))
		
		#Make the mean of the chosen curves and return
		AIF[v,:]=np.mean(chosenConccurves,0)
		AIFSI=np.mean(chosenSIcurves,0)


	return AIF

def plot_peak_baseline(dynimages, peak, slices,label):
	import numpy as np
	import roipoly
	import matplotlib.pyplot as plt

	# A function to plot out the baseline and peak signal intensities over slices
	# Extract mean baseline volume and peak volume

	peakvol=abs(dynimages[:,:,:,peak].astype('i')-dynimages[:,:,:,5].astype('i'))
	meanbase=np.zeros(len(slices))
	sdbase=np.zeros(len(slices))
	vesselpeak=np.zeros(len(slices))
	meansig=np.mean(dynimages[:,:,15,5],1)
	sdsig=np.std(dynimages[:,:,15,5],1)

	for i in range(len(slices)):
		h=plt.figure(figsize=(8,8))
		im=peakvol[slices[i],:,:,]
		#Plot the peak image and mark the vessel
		plt.imshow(im,cmap='gray',interpolation='nearest') # Show the image
		vesselpos=plt.ginput(1,timeout=100, show_clicks=True,mouse_add=1, mouse_pop=3, mouse_stop=2)
		vesselx=np.int(np.round(vesselpos[0][0]))
		vessely=np.int(np.round(vesselpos[0][1]))
		#Extract the curve for this vessel
		vesselcurve=np.squeeze(dynimages[slices[i],vessely,vesselx,:])
		#Extract maximum, mean baseline and sd in baseline
		vesselpeak[i]=np.max(vesselcurve)
		meanbase[i]=np.mean(vesselcurve[0:10])
		sdbase[i]=np.std(vesselcurve[0:10])
		plt.close()

	h,axarr=plt.subplots(3,1,sharex='all',sharey='none',figsize=(10,8),tight_layout='True')
	axarr[0].errorbar(slices,meanbase,sdbase,fmt='g.')
	axarr[0].plot(slices,vesselpeak,'x')
	axarr[0].errorbar(np.arange(len(meansig)),meansig,sdsig,fmt='b.',alpha=0.5)
	axarr[0].text(0,-20,label)
	axarr[1].autoscale='False'
	axarr[1].imshow(np.max(peakvol,2).T[50:150,:],extent=(0,192,0,30))
	axarr[2].imshow(np.max(peakvol.astype('i')-dynimages[:,:,:,5].astype('i'),1).T)
			
	return meanbase, sdbase, vesselpeak	





