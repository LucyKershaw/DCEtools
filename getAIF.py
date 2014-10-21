# Module to do AIF choosing from pixels

def getAIF(dynimages, slice, TR, flip):
	# TR should be in seconds
	import numpy as np
	import matplotlib.pyplot as plt
	import pyqtgraph as pg
	import FLASH
	print('changes?')
	# First job - find AIF peak image by finding the index of the maximum pixel value within the first 50 frames
	peakframe=np.argmax(np.amax(np.amax(dynimages[:,:,slice,0:50],0),0))

	# Show image and click the vessels
	plt.imshow(dynimages[:,:,slice,peakframe])
	vesselpos=plt.ginput(n=2,timeout=100, show_clicks=True,mouse_add=1, mouse_pop=3, mouse_stop=2) # these are not integers

	# Display a 9x9 grid of AIFs from around the points
	plt.ioff()
	g,axarr=plt.subplots(9,9,sharex=False,sharey=True,figsize=[15,10],tight_layout='True')
	

	#j is columns on the plot, i is rows

	for i in range(-4,5):
		for j in range(-4,5):
			vesselx=np.int(np.round(vesselpos[0][0]))
			vessely=np.int(np.round(vesselpos[0][1]))
			#h=g.addPlot(row=i+4,col=j+4)
			#h.plot(dynimages[vessely+i,vesselx+j,slice,:])
			axarr[i+4,j+4].plot(dynimages[vessely+i,vesselx+j,slice,:])
			axarr[i+4,j+4].set_xticklabels([])
			axarr[i+4,j+4].tick_params(labelsize=8)
			axarr[i+4,j+4].set_xbound([0,dynimages.shape[3]])
			axarr[i+4,j+4].text(0.9,0.7,str(9*(i+4)+j+4),horizontalalignment='right',transform=axarr[i+4,j+4].transAxes)
	g.draw()


	# Convert to concentration and plot again
	g,axarr=plt.subplots(9,9,sharex=False,sharey=True,figsize=[15,10],tight_layout='True')
	
	for i in range(-4,5):
		for j in range(-4,5):
			vesselx=np.int(np.round(vesselpos[0][0]))
			vessely=np.int(np.round(vesselpos[0][1]))
			axarr[i+4,j+4].plot(FLASH.SI2Conc(dynimages[vessely+i,vesselx+j,slice,:],TR,flip,1.4,10,None))
			axarr[i+4,j+4].set_xticklabels([])
			axarr[i+4,j+4].tick_params(labelsize=8)
			axarr[i+4,j+4].set_xbound([0,dynimages.shape[3]])
			axarr[i+4,j+4].set_ybound([-1,30])
			axarr[i+4,j+4].text(0.9,0.7,str(9*(i+4)+j+4),horizontalalignment='right',transform=axarr[i+4,j+4].transAxes)
	g.draw()
	plt.ion()

	# Now choose which subplots to use
	ids=input("Input the id number for the chosen AIF plots, separated by commas: ")
	ids=ids.split(',')
	ChosenAIFs=[int(item) for item in ids]
	# Work out i,js
	AIFij=[[int((x-np.remainder(x,9))/9),np.remainder(x,9)] for x in ChosenAIFs]

	chosenSIcurvesL=np.zeros((len(ChosenAIFs),dynimages.shape[3]))
	# And extract them
	for p in range(0,len(AIFij)):
		chosenSIcurvesL[p,:]=np.squeeze(dynimages[vessely+AIFij[p][0],vesselx+AIFij[p][1],slice,:])

	return chosenSIcurvesL

