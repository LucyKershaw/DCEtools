

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def makemovie(imstack):
	global x
	fig = plt.figure()
	x=0
	im = plt.imshow(imstack[:,:,x])
	#plt.ion()


	def updatefig(*args):
	    global x
	    im.set_array(imstack[:,:,x])
	    plt.title(x)
	    x+=1
	    x%=imstack.shape[2]
	    return im,

	ani = animation.FuncAnimation(fig, updatefig, interval=500)
	return ani
	#plt.show()


