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
		
		#Make the mean of the chosen curves and return
		AIF[v,:]=np.mean(chosenConccurves,0)
		AIFSI=np.mean(chosenSIcurves,0)


	return AIF

def plot_peak_baseline(dynimages, peak, slices):
	import numpy as np
	import roipoly
	import matplotlib.pyplot as plt

	# A function to plot out the baseline and peak signal intensities over slices
	# Extract mean baseline volume and peak volume

	baseline=np.mean(dynimages[:,:,:,1:10],3)
	peakvol=dynimages[:,:,:,peak]
	maxbase=np.zeros(len(slices))
	maxpeak=np.zeros(len(slices))

	for i in range(len(slices)):
		h=plt.figure()
		im=peakvol[slices[i],:,:,]
		base=baseline[slices[i],:,:]
		plt.imshow(im,cmap='gray',interpolation='nearest') # Show the image
		roi=roipoly.roipoly.roipoly() # Mark the vessel
		input('press enter to continue')
		mask=roi.getMask(im) # Make the mask

		maxpeak[i]=np.max(im*mask)# Find the peak peak value wihtin the vessel
		location=np.argmax(im*mask)
		maxbase[i]=np.ndarray.flatten(base)[location]# Find the baseline value of this pixel
		


	return maxbase, maxpeak			






