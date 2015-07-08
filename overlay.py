# Little function to overlay a mask onto a background image
# Zero pixels (where no fit was done) completely transparent and map colours set to custom transparency
# Map is interpolated to same size as background image

def overlay(background,map,maptransparency,maplimit):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	import scipy.misc
	import numpy as np

	infmask=np.isinf(map)
	nanmask=np.isnan(map)
	map[infmask]=0
	map[nanmask]=0


	#print(np.amax(map))
	#print(np.amin(map))

	#Rescale the map so that maplimit pixels are the maximum
	map[0,0]=maplimit
	clipmap=map>maplimit
	map[clipmap]=maplimit
	#Resize the map
	newmap=scipy.misc.imresize(map,(background.shape),interp='nearest',mode='L')

	#print(np.amax(newmap))
	#print(np.amin(newmap))

	#Set up the cmap for the map
	my_cmap=cm.jet
	my_cmap.set_under('k',alpha=0)
	my_cmap.set_bad('k',alpha=0)

	#Make a new figure with greyscale background
	plt.figure()
	plt.imshow(background,cmap=cm.gray,vmin=0,vmax=np.amax(background)/2,interpolation='nearest')

	#Overlay the map with transparency
	plt.imshow(newmap, cmap=my_cmap, interpolation='nearest', clim=[0.001,255],alpha=maptransparency)