# Module to select an elliptical ROI from an image stack and return a mask

from matplotlib.patches import Ellipse
from matplotlib.widgets import EllipseSelector
import matplotlib.pyplot as plt
from pylab import *
import numpy as np

def draw_ellipse(imagestack):
	#### Main function to draw ellipse on the stack
	# Make an array of point coordinates to represent the pixels in one slice
	global images, fig, ax, ellipsepos, currentslice
	images=imagestack
	mask=np.zeros_like(imagestack)
	ellipsepos=zeros((4,images.shape[2]))
	currentslice=0
	
	#Make a figure with figure and axes handles:
	fig=plt.figure(figsize=(9,9))
	ax=plt.gca()
	update_figure()
	#Start the toggle selector within these axes
	toggle_selector.ES = EllipseSelector(ax, onselect, drawtype='line', interactive='True')
	connect('key_press_event', toggle_selector)
	print('Either skip slice (press n), Draw or edit the ellipse for this slice then press d for next slice')
	input('Press any key to continue')

	print('back in draw_ellipse')
	print(ellipsepos)

	#Make a coordinate array to match one slice 
	y_int=np.arange(images.shape[0])
	x_int=np.arange(images.shape[1])
	g=np.meshgrid(x_int, y_int)
	coords=list(zip(*(c.flat for c in g)))

	#Now loop through slices, make an ellipse patch, find image coords that lie within it and set mask pixels
	for s in range(images.shape[2]):
		if sum(ellipsepos[:,s]): #if an ellipse was defined for this slice
			startx, starty, endx, endy=ellipsepos[:,s]		
			ell=Ellipse(xy=[(endx+startx)/2,(endy+starty)/2],width=endx-startx, height=endy-starty, angle=0)
			ellipsepoints=np.vstack([p for p in coords if ell.contains_point(p,radius=0)])
			mask[ellipsepoints[:,1],ellipsepoints[:,0],s]=1

	ellipsepoints=0
	pos=0
	fig=0
	ax=0
	return mask

def onselect(eclick, erelease): # define function to update pos with the start and end click positions
    #eclick and erelease are matplotlib events at press and release'
    global pos
    pos=np.array([eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata])
    print(' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata))
    print(' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata))
    	   		

def toggle_selector(event): # define function for key press
	global ellipsepos, currentslice, editing
	if event.key in ['D', 'd'] and toggle_selector.ES.active:
		#When next slice is requested, save the current startpos and endpos in the ellipsepos variable
		ellipsepos[:,currentslice]=pos
		# Increment the current slice and.....
		currentslice=currentslice+1
		if currentslice>images.shape[2]-1: #If this takes us beyond the last slice
			plt.close(fig) #close the figure and return
			return
		else: #otherwise, update the figure
			update_figure()
	if event.key in ['N','n'] and toggle_selector.ES.active:
		#If press s for skip slice, just increment the currentslice and update the figure or return
		currentslice=currentslice+1
		if currentslice>images.shape[2]-1:
			plt.close(fig)
			return
		else:
			update_figure()



def update_figure():
		ax.imshow(images[:,:,currentslice])
		fig.suptitle(str(currentslice))
		fig.canvas.draw()

























