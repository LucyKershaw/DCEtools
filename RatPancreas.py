# Code to analyse rat pancreas data from T1DM study

#Imports
import numpy as np
import dicom
import matplotlib.pyplot as plt
import os
import glob
import shutil
from tkinter import Tk
from tkinter import filedialog
from roipoly import roipoly
from scipy.optimize import curve_fit

import scipy.misc
import scipy.optimize
import pprint


class ratdata(object):
	def __init__(self,ratdirect):
		datadirect=os.path.join(ratdirect,'data')
		self.datadirect=datadirect
		# Make a directory to hold the analysis results
		# Add an Analysis folder in the main directory, if it doesn't already exist
		if not os.path.isdir(os.path.join(os.path.split(datadirect)[0],'Analysis')):
			print('Making analysis directory')
			os.mkdir(os.path.join(os.path.split(datadirect)[0],'Analysis'))

	def disp_seriesfolders(self):
		subdirects=[x for x in os.listdir(self.datadirect) if os.path.isdir(os.path.join(self.datadirect, x))]
		for item in subdirects:
			print(item)

	def read_structural(self,seriestag): #method to read image stack - read in by file name due to poor dicom conversion
		# Find folder that matches seriestag
		stackfolder=glob.glob(os.path.join(self.datadirect,seriestag))
		if not stackfolder:
			print('Folder not found')
			return
		print('Found ', stackfolder[0])

		# find .dcm files
		filenames=glob.glob(os.path.join(stackfolder[0],'*.dcm'))
		print('Check file ordering: ')
		for item in filenames:
			print(os.path.split(item)[1])
		
		numfiles=len(filenames)
		print("Reading "+str(numfiles)+" files")
		
		# read the first file to find out size, and add pixel size info to T2winfo structure
		info=dicom.read_file(filenames[0])
		self.structuralinfo=np.zeros(1,dtype=[('pixelsize','f8')])
		self.structuralinfo['pixelsize']=float(info.PixelSpacing[0])
		im=info.pixel_array

		ims=np.zeros(np.array([im.shape[0],im.shape[1],numfiles]))
		for i in range(0,numfiles):
			temp=dicom.read_file(filenames[i])
			ims[:,:,i]=np.flipud(temp.pixel_array)
		self.structural=ims

	def read_dynamics(self,seriestag): #Method to read in single slice dynamic images, by filename due to poor dicom conversion
		#Look for dynamics.npy in the Analysis folder first
		if os.path.isfile(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dynamics.npy')):
			print('reading from saved array')
			self.dynims=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dynamics.npy'))
			self.dyninfo=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dyninfo.npy'))
			return
		# Use known series tag by calling disp_seriesfolders first
		
		# Find folder that matches seriestag
		DynFolder=glob.glob(os.path.join(self.datadirect,seriestag))
		if not DynFolder:
			print('Folder not found')
			return
		print('Found ',DynFolder[0])

		# Find all the dynamic filenames
		dynfiles=glob.glob(os.path.join(self.datadirect,DynFolder[0],'*.dcm'))
		numfiles=len(dynfiles)
		print("Reading "+str(numfiles)+" dynamic files")

		# Read the last file to work out size, check for manufacturer and fill out info variables
		info=dicom.read_file(dynfiles[-1])
		im=info.pixel_array
		# Find pixel sizes and other required dynamic info
		self.dyninfo=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8'),('tres','f8'),('numtimepoints','i4'),('numslices','i4')])
		self.dyninfo['pixelsize']=float(info.PixelSpacing[0])
		self.dyninfo['TR']=float(info.RepetitionTime)
		self.dyninfo['FlipAngle']=float(info.FlipAngle)

		# Make an array to hold the dynamic data
		dynims=np.zeros((im.shape[0],im.shape[1],numfiles),dtype='uint16')
		
		# Read files into the array
		for i in range(0,len(dynfiles)):
			temp=dicom.read_file(dynfiles[i]) # Read file
			dynims[:,:,i]=np.flipud(temp.pixel_array) # Read into the right part of the array

		# save this file as an npy array for next time
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dynamics.npy'),dynims)
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dyninfo.npy'),self.dyninfo)
		print("dynamic image array size is "+str(dynims.shape))
		self.dynims=dynims

	def read_MOLLI(self,seriestag,PreOrPost): #Method to read in single slice MOLLI images, by filename due to poor dicom conversion
		#Look for MOLLIpre and MOLLIpost.npy in the Analysis folder first
		if PreOrPost=='Pre' and os.path.isfile(os.path.join(os.path.split(self.datadirect)[0],'Analysis','MOLLIpre.npy')):
			print('reading from saved array')
			self.MOLLIpre=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','MOLLIpre.npy'))
			#self.dyninfo=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dyninfo.npy'))
			return

		if PreOrPost=='Post' and os.path.isfile(os.path.join(os.path.split(self.datadirect)[0],'Analysis','MOLLIpost.npy')):
			print('reading from saved array')
			self.MOLLIpost=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','MOLLIpost.npy'))
			#self.dyninfo=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dyninfo.npy'))
			


		# Use known series tag by calling disp_seriesfolders first
		# Find folder that matches seriestag
		MOLLIFolder=glob.glob(os.path.join(self.datadirect,seriestag))
		if not MOLLIFolder:
			print('Folder not found')
			return
		print('Found ',MOLLIFolder[0])

		# Find all the filenames
		MOLLIfiles=glob.glob(os.path.join(self.datadirect,MOLLIFolder[0],'*.dcm'))
		numfiles=len(MOLLIfiles)
		print("Reading "+str(numfiles)+" MOLLI files")

		# Read the last file to work out size, check for manufacturer and fill out info variables
		info=dicom.read_file(MOLLIfiles[-1])
		im=info.pixel_array
		# Find pixel sizes and other required dynamic info
		# self.dyninfo=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8'),('tres','f8'),('numtimepoints','i4'),('numslices','i4')])
		# self.dyninfo['pixelsize']=float(info.PixelSpacing[0])
		# self.dyninfo['TR']=float(info.RepetitionTime)
		# self.dyninfo['FlipAngle']=float(info.FlipAngle)

		# Make an array to hold the dynamic data
		MOLLIims=np.zeros((im.shape[0],im.shape[1],numfiles),dtype='uint16')
		
		# Read files into the array
		for i in range(0,len(MOLLIfiles)):
			temp=dicom.read_file(MOLLIfiles[i]) # Read file
			MOLLIims[:,:,i]=np.flipud(temp.pixel_array) # Read into the right part of the array

		# save this file as an npy array for next time

		if PreOrPost=='Pre':
			np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','MOLLIpre.npy'),MOLLIims)
			#np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dyninfo.npy'),self.dyninfo)
			self.MOLLIpre=MOLLIims

		if PreOrPost=='Post':
			np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','MOLLIpost.npy'),MOLLIims)
			#np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dyninfo.npy'),self.dyninfo)
			self.MOLLIpost=MOLLIims

	def calc_percentenh(self):
		#Method to calculate percentage enhancement over baseline
		baseline=np.mean(self.dynims[:,:,:5],2)
		numims=self.dynims.shape[2]
		movingav=np.zeros(self.dynims.shape)
		for i in range(2,numims-1):
			movingav[:,:,i]=np.mean(self.dynims[:,:,i-2:i+3],2)

		percentenh=np.zeros(self.dynims.shape)
		for i in range(numims):
			percentenh[:,:,i]=100*((movingav[:,:,i]-baseline)/baseline)

		self.percentenh=percentenh
		self.postcontrast=np.mean(self.dynims[:,:,-10:],2)
		self.baseline=baseline
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','percentenh.npy'),percentenh)
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','baseline.npy'),baseline)
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','postcontrast.npy'),np.mean(self.dynims[:,:,-10:],2))
		
	def calc_maxenh(self,use_saved_roi=1):
		#Method to calculate maximum enhancement for a region by finding the frame with the maximum mean over the region

		#Load required images if not already present - postcontrast and percentage enhancement
		if not hasattr(self,'percentenh'):
			print('Loading from saved percentage enhancement')
			self.percentenh=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','percentenh.npy'))

		if not hasattr(self,'postcontrast'):
			print('Loading from saved postcontrast')
			self.postcontrast=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','postcontrast.npy'))

		#Draw region on post-contrast images if needed
		if use_saved_roi==0:
			postcontrast=self.postcontrast
			plt.imshow(postcontrast)
			roi=roipoly.roipoly()
			input('press enter to continue')
			mask=roi.getMask(postcontrast)
			# Save the mask just in case
			np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','mask.npy'),mask)
			numpix=np.sum(mask)
			curve=np.mean(self.percentenh[mask==1,:],0)
			sdcurve=np.std(self.percentenh[mask==1,:],0)

		if use_saved_roi==1:
			mask=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','mask.npy'))
			numpix=np.sum(mask)
			curve=np.mean(self.percentenh[mask==1,:],0)
			sdcurve=np.std(self.percentenh[mask==1,:],0)

		if use_saved_roi==2:
			mask=self.mask
			curve=np.zeros(self.percentenh.shape[2])
			sdcurve=curve
			for i in range(0,self.percentenh.shape[2]):
				curve[i]=np.mean(self.percentenh[mask[:,:,i]==1,i])
				sdcurve[i]=np.std(self.percentenh[mask[:,:,i]==1,i])

		plt.plot(curve[2:len(curve)-1])
		#plt.plot(curve[2:len(curve)-1]-sdcurve[2:len(curve)-1])
		#plt.plot(curve[2:len(curve)-1]+sdcurve[2:len(curve)-1])
		ax=plt.gca()
		ax.set_ybound([-10,180])

		maximarg=np.argmax(curve)
		roimax=np.max(curve)
		roisd=sdcurve[maximarg]

		plt.text(75,0,'maximum {0:.2f} +- {1:.2f} at {2}'.format(roimax,roisd,maximarg))
		print('Maximum occurs at '+str(maximarg))
		print('Maximum mean percentage enhancement in roi is '+str(roimax))

		maxim=self.percentenh[:,:,maximarg]
		#plt.figure()
		#plt.imshow(maxim)
		self.maxim=maxim
		self.roimax=roimax
		self.mask=mask
		self.roisd=roisd

	def calc_sigmoid_fits(self):
		# Method to use previously saved roi with a new liver roi to normalise the enhancement curves

		#Load required images and pancreas mask if not already present
		if not hasattr(self,'dynims'):
			print('Loading dynims from file')
			self.dynims=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','dynamics.npy'))

		if not hasattr(self,'mask'):
			print('Loading from saved mask')
			self.mask=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','mask.npy'))

		if not hasattr(self,'baseline'):
			print('Loading from saved baseline image')
			self.baseline=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','baseline.npy'))

		if not hasattr(self,'postcontrast'):
			print('Loading from saved baseline image')
			self.postcontrast=np.load(os.path.join(os.path.split(self.datadirect)[0],'Analysis','postcontrast.npy'))


		#Draw region on liver
		postcontrast=self.postcontrast
		plt.imshow(postcontrast)
		roi=roipoly.roipoly()
		input('press enter to continue')
		livermask=roi.getMask(postcontrast)
		# Save the livermask just in case
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','livermask.npy'),livermask)


		#Subtract baseline voxelwise
		numims=self.dynims.shape[2]
		enhoverbase=np.zeros(self.dynims.shape)
		for i in range(numims):
			enhoverbase[:,:,i]=100*((self.dynims[:,:,i]-self.baseline)/self.baseline)

		#Calculate mean curves for liver and pancreas regions
		roicurve=np.mean(enhoverbase[self.mask==1,:],0)
		livercurve=np.mean(enhoverbase[livermask==1,:],0)
		times=np.arange(0,numims*25/60,25/60)

		#Fit to sigmoid
		def sigmoid(t, A,b,t0):
			return A/(1+np.exp(-1*b*(t-t0)))

		liveropt, livercov = curve_fit(sigmoid, times, livercurve, p0=[100,0.1,25/60])
		livererr = np.sqrt(np.diag(livercov))
		pancopt, panccov = curve_fit(sigmoid, times, roicurve, p0=[100,0.1,25/60])
		pancerr = np.sqrt(np.diag(panccov))
		normopt,normcov = curve_fit(sigmoid, times, roicurve/liveropt[0], p0=[100,0.1,25/60])
		normerr= np.sqrt(np.diag(normcov))

		plt.plot(times,roicurve,'g.')
		plt.plot(times,livercurve,'r.')
		plt.plot(times, sigmoid(times, *pancopt), 'g-')
		plt.plot(times, sigmoid(times, *liveropt), 'r-')
		
		print(pancopt[0],pancerr[0],pancopt[1],pancerr[1],liveropt[0],livererr[0],liveropt[1],livererr[1],normopt[0],normerr[0],normopt[1],normerr[1])

		self.roicurve=roicurve
		self.livercurve=livercurve

		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','roicurve.npy'),roicurve)
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','livercurve.npy'),livercurve)
		np.save(os.path.join(os.path.split(self.datadirect)[0],'Analysis','roinormcurve.npy'),roicurve/liveropt[0])

		#plt.text(75,0,'maximum {0:.2f} +- {1:.2f} at {2}'.format(roimax,Normcurvesd,maximarg))
		#print('Maximum occurs at '+str(maximarg))
		#print('Maximum mean percentage enhancement in roi is '+str(roimax))

















