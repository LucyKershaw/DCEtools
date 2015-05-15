# Class to extract contours from a structure set file
#Imports
import dicom
import numpy as np
from tkinter import Tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import glob
import os

class contourdisplay(object):
	def __init__(self,contourdirect='select'):
		self.contourdirectory=contourdirect
		if contourdirect=='select':
				# Choose directory containing contour and images if not already stated
				root = Tk()
				#root.withdraw()
				direct = filedialog.askdirectory(title="select contour and image directory")
				root.destroy()
				print(direct)
				self.contourdirectory=direct

	def get_contours(self):
		#Look for structure set file in the contour directory
		ssfilename=glob.glob(os.path.join(self.contourdirectory,'*RS*'))
		if len(ssfilename)>1:
			print('More than one structure set, choose the correct one')
			root = Tk()
			ssfilename = filedialog.askopenfilename(title="select structure set file")
			root.destroy()
		print(ssfilename)
		self.ssfilename=ssfilename[0]
		#Read the structure set file
		self.ss=dicom.read_file(ssfilename[0])

		#Find number of contours in structure set
		self.numcontours=len(self.ss.ROIContours)

		#Get contour names and associated numbers
		self.ROINames=[self.ss.StructureSetROIs[x].ROIName for x in range(0,self.numcontours)]
		self.ROINums=[self.ss.StructureSetROIs[y].ROINumber for y in range(0,self.numcontours)]
		self.ROIContournums=[self.ss.ROIContours[z].ReferencedROINumber for z in range(0,self.numcontours)]

	def get_images(self):
		#Read images in from the same folder
		imfilenames=glob.glob(os.path.join(self.contourdirectory,'*MR*'))
		IPP=[dicom.read_file(x).ImagePositionPatient for x in imfilenames]
		IPPz=[f[2] for f in IPP]
		together=zip(IPPz,imfilenames,IPP)
		together=sorted(together)
		sortedimfilenames=[y[1] for y in together]
		sortedUIDfromimfilenames=[x.split('MR.')[1] for x in sortedimfilenames]

		sortedIPP=[z[2] for z in together]
		sortedIPPz=[p[0] for p in together]

		#Read in the files in the proper order
		#Read in first one to get information needed
		info=dicom.read_file(sortedimfilenames[0])
		cols=info.Columns
		rows=info.Rows
		PS=info.PixelSpacing
		im=info.pixel_array

		images=np.zeros(np.array([im.shape[0],im.shape[1],len(sortedimfilenames)]))
		
		for i in range(0,len(sortedimfilenames)):
			temp=dicom.read_file(sortedimfilenames[i])
			images[:,:,i]=temp.pixel_array
		
		self.images=images
		self.sortedIPP=sortedIPP
		self.sortedUIDfromimfilenames=sortedUIDfromimfilenames
		self.rows=rows
		self.cols=cols
		self.PS=PS

	def disp_contour(self):
		#Method to display one or more contours on associated images
		#Choose contours for display

		[print(item) for item in enumerate(self.ROINames)]
		ids=input("Input the id number for the chosen contour(s), separated by commas: ")
		ids=ids.split(',')
		ROIchoice=[self.ROINames[int(item)] for item in ids]
		print("Chosen rois:")
		[print(item) for item in ROIchoice]


		#Check that ROIContournum matches ROINum
		for n in range(0,len(ids)):
			if self.ROINums[int(ids[n])]!=self.ROIContournums[int(ids[n])]:
				print(self.ROINames[int(ids[n])])
				print('Warning - contours stored out of order, check this is the right structure')

		#For each slice, find any matching contour slices and plot
		plt.ioff()
		for i in range(0,np.shape(self.images)[2]):
			#Find slice UID and display image
			imUID=self.sortedUIDfromimfilenames[i]
			currentfig=plt.figure(str(i))
			#Plot the image with z order of 1
			plt.imshow(self.images[:,:,i],extent=[self.sortedIPP[i][0],self.sortedIPP[i][0]+self.cols*self.PS[0],self.sortedIPP[i][1]+self.rows*self.PS[0],self.sortedIPP[i][1]],interpolation='nearest',cmap='gray',zorder=1)

			numplots=0
			#For each contour selected, find if this UID appears in the list
			linecols=['r','b','g','y','c','k','w','m']
			for m in range(0,len(ids)):
				Contour=self.ss.ROIContours[int(ids[m])].ContourSequence
				ContourRefdUID=[Contour[x].ContourImageSequence[0].ReferencedSOPInstanceUID for x in range(0,len(Contour))] #Find all image UIDs referred to
				if imUID in ContourRefdUID: #if this image appears in the contour refs
					ind=ContourRefdUID.index(imUID)
					#plot contour using correct coordinates, and with z order to leave space for image shifts
					SliceContourData=Contour[ind].ContourData
					plt.plot(SliceContourData[0::3],SliceContourData[1::3],linecols[m],label=ROIchoice[m],zorder=100+m)
					numplots=numplots+1
			# #If this slice has no contours on it, close the figure
			if numplots==0: 
				plt.close(currentfig)
			else:
				plt.legend()
				currentfig.show()

		plt.ion()

	def shift_contour(self,distance):
		#Find numbers for current figures
		openfiglabels=plt.get_figlabels()
		openfignums=plt.get_fignums()
		#for each one, show the image again but shift in the AP direction and replot between the original image and the lines
		for i in range(len(openfiglabels)):
			plt.figure(openfignums[i])
			plt.imshow(self.images[:,:,openfiglabels[int(i)]],extent=[self.sortedIPP[i][0],self.sortedIPP[i][0]+self.cols*self.PS[0],self.sortedIPP[i][1]+self.rows*self.PS[0]+distance,self.sortedIPP[i][1]+distance],interpolation='nearest',cmap='gray',zorder=2)			
