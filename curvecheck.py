# Helper function to click on a pixel and plot the curve for that pixel

import matplotlib.pyplot as plt
import numpy as np

def curvecheck(imstack):
	h=plt.figure()
	plt.imshow(imstack[:,:,0])
	pix_pos=plt.ginput(n=1,timeout=100, show_clicks=True,mouse_add=1, mouse_pop=3, mouse_stop=2) # these are not integers
	print(pix_pos)
	pix_x=np.int(np.round(pix_pos[0][0]))
	pix_y=np.int(np.round(pix_pos[0][1]))

	j=plt.figure()
	plt.plot(np.squeeze(imstack[pix_y,pix_x,:]))

	return np.squeeze(imstack[pix_y,pix_x,:])