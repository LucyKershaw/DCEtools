# Main file for defining the patient class to do DCE processing pixel by pixel

# Imports
import numpy as np
import dicom
import matplotlib.pyplot as plt
import os
import glob
import shutil
import csv
from tkinter import Tk
from tkinter import filedialog
from roipoly import roipoly

import scipy.misc
import scipy.signal
import SRTurboFLASH
import IRTurboFLASH
import IRTrueFISP
import scipy.optimize
import scipy.ndimage.morphology
import TwoCXM
import TwoCUM
import AATH
import ExtKety
import FLASH
import getAIF
import pprint
import overlay

class patient(object): # patient inherits from the object class
	
	def __init__(self,patientdirectory='select'):
		self.patientdirect=patientdirectory
		self.dicomdirect=os.path.join(self.patientdirect,'DICOM')
		if patientdirectory=='select':
			# Choose patient directory if not already stated
			root = Tk()
			#root.withdraw()
			direct = filedialog.askdirectory(initialdir='C:/Data/',title="select patient directory")
			root.destroy()
			print(direct)
			self.patientdirect=direct
			direct=os.path.join(direct,'DICOM')
			self.dicomdirect=direct

		if os.path.isfile(os.path.join(self.patientdirect,'Analysis','hct.npy')):
			self.hct=np.load(os.path.join(self.patientdirect,'Analysis','hct.npy'))
		else: 
			print('No hct found in patient folder, using default value of 0.4')
			self.hct=0.4

	# method for displaying data attributed
	def disp_attr(self):
		print("Images:")
		if hasattr(self, 'T2wims'):
			print("    T2w images: "+str(self.T2wims.shape))
			print("    "+str(self.T2winfo.dtype.names))
			print("    "+str(self.T2winfo))
			print("")
		if hasattr(self, 'dynims'):
			print("    Dynamic images: "+str(self.dynims.shape))
			print("    "+str(self.dyninfo.dtype.names))
			print("    "+str(self.dyninfo))
			print("")
		if hasattr(self,'T1data'):
			print("     T1 measurement images: "+str(self.T1data.shape))
			print("    "+str(self.T1info.dtype.names))
			print("    "+str(self.T1info))

	# Method for showing available series folders - use to find suitable tag when reading images
	def disp_seriesfolders(self):
		subdirects=[x for x in os.listdir(self.dicomdirect) if os.path.isdir(os.path.join(self.dicomdirect, x))]
		for item in subdirects:
			print(item)
		#pprint.pprint(subdirects)

	# Define methods to read in all the required images
	#####################################################
	def read_T2w(self,seriestag): #method to read the T2w image files
		# Use known series tag by calling disp_seriesfolders first
		
		# Find folder that matches seriestag
		T2wFolder=glob.glob(os.path.join(self.dicomdirect,seriestag))
		if not T2wFolder:
			print('Folder not found')
			return
		print('Found ',T2wFolder[0])

		# find .dcm files
		filenames=glob.glob(os.path.join(T2wFolder[0],'*.dcm'))
		#print(filenames)
		numfiles=len(filenames)
		print("Reading "+str(numfiles)+" files")
		# read the first file to find out size, and add pixel size info to T2winfo structure
		info=dicom.read_file(filenames[0])
		self.T2winfo=np.zeros(1,dtype=[('pixelsize','f8')])
		self.T2winfo['pixelsize']=float(info.PixelSpacing[0])
		im=info.pixel_array

		T2wims=np.zeros(np.array([im.shape[0],im.shape[1],numfiles]))
		for i in range(0,numfiles):
			temp=dicom.read_file(filenames[i])
			imnum=temp.InstanceNumber
			T2wims[:,:,imnum-1]=temp.pixel_array
		self.T2wims=T2wims


	def read_dynamics(self,seriestag):
		
		#Look for dynamics.npy in the Analysis folder first
		if os.path.isfile(os.path.join(self.patientdirect,'Analysis','dynamics.npy')):
			print('reading from saved array')
			self.dynims=np.load(os.path.join(self.patientdirect,'Analysis','dynamics.npy'))
			self.dyninfo=np.load(os.path.join(self.patientdirect,'Analysis','dyninfo.npy'))
			self.t=np.arange(0,self.dyninfo['tres'][0]*self.dyninfo['numtimepoints'][0],self.dyninfo['tres'][0])

			return
		# Use known series tag by calling disp_seriesfolders first
		
		# Find folder that matches seriestag
		DynFolder=glob.glob(os.path.join(self.dicomdirect,seriestag))
		if not DynFolder:
			print('Folder not found')
			return
		print('Found ',DynFolder[0])

		# Find all the dynamic filenames
		dynfiles=glob.glob(os.path.join(self.dicomdirect,DynFolder[0],'*.dcm'))
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

		# Detect Siemens dicom acquired as separate series or Siemens dicom acquired as one series
		if hasattr(info,"TemporalPositionIdentifier"): # If this tag exists...
			print(info.TemporalPositionIdentifier)
			if info.TemporalPositionIdentifier != '': # And it's not empty, the job is easy...
				print('TemporalPositionIdentifier found')
				numtimepoints=int(info.TemporalPositionIdentifier)
				self.dyninfo['numtimepoints']=numtimepoints
				self.dyninfo['tres']=round(float(info.AcquisitionDuration/numtimepoints),3)
				numslices=int(len(dynfiles)/numtimepoints)
				self.dyninfo['numslices']=numslices
			else: #If the tag exists but it's empty, detect the number of series from the names
				print('TemporalPositionIdentifier empty, using individual series numbers')
				filenamesonly=[os.path.basename(Y) for Y in dynfiles]
				dynseriesnums=set([f.split('_')[0] for f in filenamesonly])
				numtimepoints=len(dynseriesnums)
				self.dyninfo['numtimepoints']=numtimepoints
				numslices=np.int(np.ceil(len(dynfiles)/numtimepoints))
				self.dyninfo['numslices']=numslices			
				print('WARNING, set tres manually - cannot calculate properly from headers')

		else: # If the tag doesn't exist, get numtimepoints from AcquisitionNumber of last file
			print('No temporal position identifier tag, using AcquisitionNumber')
			numtimepoints=int(info.AcquisitionNumber)
			self.dyninfo['numtimepoints']=numtimepoints
			try:
				self.dyninfo['tres']=float(info[0x0051,0x100a].value.split('TA ')[1])
			except TypeError:
				self.dyninfo['tres']=float(info[0x0051,0x100a].value.decode('ASCII').split('TA ')[1]) # Added to take account of new anon using OB rather than LO
			numslices=int(len(dynfiles)/numtimepoints)
			self.dyninfo['numslices']=numslices

		# Make an array to hold the dynamic data
		dynims=np.zeros(np.array([im.shape[0],im.shape[1],numslices,numtimepoints]),dtype='uint16')
		
		# Read files into the array
		for i in range(0,len(dynfiles)):
			temp=dicom.read_file(dynfiles[i]) # Read file
			if hasattr(temp,"TemporalPositionIdentifier"): #if this exists,
				if temp.TemporalPositionIdentifier!='': #and it isn't empty
					temporalpos=temp.TemporalPositionIdentifier #just use it
				else:
					# but if not, use the AcquisitionNumber instead
					temporalpos=temp.AcquisitionNumber
			else:
				temporalpos=temp.AcquisitionNumber #And if it doesn't have TemporalPositionIdentifier, use the Acquisition Number
			timept=temporalpos-1 # Find temporal position (begin at zero!)
			
			if temp.Manufacturer=='Philips Medical Systems':
				slice=int(np.floor(i/numtimepoints)) # Work out slice number for Philips (begin at zero!) 
				dynims[:,:,slice,timept]=temp.pixel_array # Read into the right part of the array
				#print(i,slice, timept, temp.SliceLocation, '0020, 0012', temp[0x0020,0x0012].value, 'Inst num, 0020, 0013',temp[0x0020,0x0013].value)
			if temp.Manufacturer=='SIEMENS':
				slice=int(i-timept*numslices) # Work out the slice number for Siemens
				#print(i,slice, timept, temp.SliceLocation, '0020, 0012', temp[0x0020,0x0012].value, 'Inst num, 0020, 0013',temp[0x0020,0x0013].value)
				dynims[:,:,slice,timept]=temp.pixel_array # Read into the right part of the array

		# save this file as an npy array for next time
		np.save(os.path.join(self.patientdirect,'Analysis','dynamics.npy'),dynims)
		np.save(os.path.join(self.patientdirect,'Analysis','dyninfo.npy'),self.dyninfo)
		print("dynamic image array size is "+str(dynims.shape))
		self.t=np.arange(0,self.dyninfo['tres'][0]*self.dyninfo['numtimepoints'][0],self.dyninfo['tres'][0])
		self.dynims=dynims

	def read_T1data(self, seriestag, measnum=0, usefolders=range(0,6,1), VFA=0, multimeas=1):
	 #Set VFA flag to 1 if necessary, multimeas is for multiple dynamics in a VFA
		T1Folders=glob.glob(os.path.join(self.dicomdirect,seriestag))
		if not T1Folders:
			print('Folder not found')
			return
		print('Found:')
		pprint.pprint(T1Folders)
		T1Folders=[T1Folders[x] for x in usefolders]
		print('Using these folders:' )
		pprint.pprint(T1Folders)
		if measnum==0:
			print('This is going to be saved as the primary, precontrast T1 measurement')
		if measnum==1:
			print('This is going to be treated as a secondary T1 measurement')


		# in the first directory, read in one file to check sizes etc
		T1files=glob.glob(os.path.join(T1Folders[0],'*.dcm')) #T1files is number of files in the series
		info=dicom.read_file(T1files[0])
		# print(T1files[0])
		im=info.pixel_array		
		# plt.imshow(im)
		T1info=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8',len(T1Folders)),('TIs','f8',len(T1Folders)),('N','f8')])
		T1info['pixelsize']=float(info.PixelSpacing[0])
		T1info['TR']=float(info.RepetitionTime)
		T1info['N']=int(info.EchoTrainLength)
		T1info['FlipAngle'][0][0]=float(info.FlipAngle)
		self.T1info=T1info

		if info.Manufacturer=='SIEMENS':
			print('Warning - values for N and TR may be incorrect for Siemens dicom - please check:')
			print('N = ' + str(self.T1info['N']))
			print('TR = '+str(self.T1info['TR']))
			# self.T1info['N']=0
			# self.T1info['TR']=0


		# Make an array of the right size to hold the images
		T1ims=np.zeros(np.array([im.shape[0],im.shape[1],len(T1files),len(T1Folders)]))
		
		# Read in the rest, checking TI or flip angle for each folder (don't rely on series names if Philips)
		for y in range(0,len(T1Folders)): # for each subfolder
			T1files=glob.glob(os.path.join(self.dicomdirect,T1Folders[y],'*.dcm')) # find dicom files
			temp=dicom.read_file(T1files[0]) # read first one
			if VFA==0:
				T1info['TIs'][0][y]=temp[0x0018,0x0082].value # set the TI value in T1info
			if VFA==1:
				T1info['FlipAngle'][0][y]=float(temp.FlipAngle)
			
			ss=1
			si=0
			if info.Manufacturer=='Philips Medical Systems':
				#If Philips, make sure the scaling of the images is done
				ss=float(temp[0x2005,0x100e].value) #extract the scale slope
				si=float(temp[0x2005,0x100d].value) #extract the scale intercept

			T1ims[:,:,0,y]=(temp.pixel_array-si)/ss

			for z in range(1,len(T1files)): # for the rest of the files
				temp=dicom.read_file(T1files[z])
				ss=1
				si=0
				if temp.Manufacturer=='Philips Medical Systems':
					ss=float(info[0x2005,0x100e].value) #extract the scale slope
					si=float(info[0x2005,0x100d].value) #extract the scale intercept

				T1ims[:,:,z,y]=(temp.pixel_array-si)/ss

		# Make an average image for the number of multiple measurements specified (usually 1)
		# For Siemens, order is {slice 0, slice 1... slice n}{slice 0, slice 1... slice n}
		
		if multimeas!=1:
			if temp.Manufacturer=='Philips Medical Systems':
				print('Multiple measurements for Philips not yet implemented')
				return
			T1imsAv=np.zeros(np.array([im.shape[0],im.shape[1],int(len(T1files)/multimeas),len(T1Folders)]))
			numslices=int(len(T1files)/multimeas)
			for k in range(numslices):
				T1imsAv[:,:,k,:]=np.mean(T1ims[:,:,k+numslices::numslices,:],2)

			T1ims=T1imsAv

		print("T1 measurement image array size is "+str(T1ims.shape))
		if VFA==0:
			print('Inversion times are '+str(T1info['TIs']))
		if VFA==1:
			print('Flip angles are '+str(T1info['FlipAngle']))
		if measnum==0:
			self.T1data=T1ims
			self.T1info=T1info
		if measnum==1:
			self.T1data2=T1ims
			self.T1info2=T1info

	# Finally, a method to read all the images
	def read_ims(self,T2wseriestag,dynseriestag,T1seriestag):
		self.read_T2w(T2wseriestag)
		self.read_T1data(T1seriestag)
		self.read_dynamics(dynseriestag)

	# methods for AIF
	#####################################################

	def read_AIF_fromfile(self):
		# read existing AIF from file in patient/Analysis directory
		print('looking for AIF file...')
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','AIF.npy')):
			print('AIF file not found')
			return	
		print('Found')
		AIF=np.load(os.path.join(self.patientdirect,'Analysis','AIF.npy'))
		self.AIF=np.squeeze(AIF)

	def get_pixel_AIF(self, usesag=0, numvessels=2, slicegroup=[3,4,5,6,7,8,9,10,11],blinepts=15):
		# Display dynamics and choose slices from which AIF should be extracted
		# First job - find AIF peak image by finding the index of the maximum pixel value within the first 50 frames for the central slice
		dynims=np.copy(self.dynims)
		midslice=np.int(np.floor(dynims.shape[2]/2))
		#slicegroup=[3,4,5,6,7,8,9,10,11]
		peakframe=np.argmax(np.amax(np.amax(dynims[:,:,midslice,0:50],0),0))
		print(peakframe)
		if usesag==1:
			midslice=np.int(np.floor(dynims.shape[1]/2))
			peakframe=21
			#peakframe=np.argmax(np.max(dynims[:,midslice,:,0:50],(0,1)))	
			print(peakframe)		
			dynims=np.swapaxes(np.swapaxes(dynims,0,2),0,1)


		h2=plt.figure()
		plt.subplot(3,3,1)
		plt.imshow(dynims[:,:,slicegroup[0],peakframe])
		plt.title(slicegroup[0])
		plt.subplot(3,3,2)
		plt.imshow(dynims[:,:,slicegroup[1],peakframe])
		plt.title(slicegroup[1])
		plt.subplot(3,3,3)
		plt.imshow(dynims[:,:,slicegroup[2],peakframe])
		plt.title(slicegroup[2])
		plt.subplot(3,3,4)
		plt.imshow(dynims[:,:,slicegroup[3],peakframe])
		plt.title(slicegroup[3])
		plt.subplot(3,3,5)
		plt.imshow(dynims[:,:,slicegroup[4],peakframe])
		plt.title(slicegroup[4])
		plt.subplot(3,3,6)
		plt.imshow(dynims[:,:,slicegroup[5],peakframe])
		plt.title(slicegroup[5])
		plt.subplot(3,3,7)
		plt.imshow(dynims[:,:,slicegroup[6],peakframe])
		plt.title(slicegroup[6])
		plt.subplot(3,3,8)
		plt.imshow(dynims[:,:,slicegroup[7],peakframe])
		plt.title(slicegroup[7])
		plt.subplot(3,3,9)
		plt.imshow(dynims[:,:,slicegroup[8],peakframe])
		plt.title(slicegroup[8])



		ids=input("Input the id number for the slices to be used, separated by commas: ")
		ids=ids.split(',')
		ChosenSlices=[int(item) for item in ids]
		plt.close(h2)

		# Get resulting AIFs for these slices
		AIFtemp=np.zeros((numvessels,dynims.shape[3],len(ChosenSlices)))
		for i in range(0,len(ChosenSlices)):
			AIFtemp[:,:,i]=getAIF.getAIF(dynims,ChosenSlices[i], self.dyninfo['TR']/1000, self.dyninfo['FlipAngle'],peakframe,numvessels,blinepts)

		# Choose from plot which one to use
		h3=plt.figure()
		for ii in range(0,len(ChosenSlices)): 
			plt.plot(np.squeeze(AIFtemp[0,:,ii]),label=str(ii)+' - '+str(ChosenSlices[ii])+'L')
			if numvessels==2:
				plt.plot(np.squeeze(AIFtemp[1,:,ii]),label=str(ii)+' - '+str(ChosenSlices[ii])+'R')
			plt.legend()
		
		ChooseAIF=input('Which AIF should be used? Input index number, comma, 0 for L, 1 for R OR Av for average all')
		if ChooseAIF=='Av':
			AIF=np.mean(AIFtemp,2)
		else:
			ChooseAIF=ChooseAIF.split(',')
			AIF=AIFtemp[int(ChooseAIF[1]),:,int(ChooseAIF[0])]
		np.save(os.path.join(self.patientdirect,'Analysis','AIFall.npy'),AIFtemp)
		np.save(os.path.join(self.patientdirect,'Analysis','AIF.npy'),AIF) #save chosen AIF in Analysis folder as AIF.
		self.AIF=AIF
		
	# Initial processing
	#####################################################
	def make_mask(self, use_dyn=0, tightmask=0, use_roipoly=1, save=1):
		#Method to make a binary mask covering the area of interest - use T2w images
		
		if use_dyn==1:
			if not hasattr(self,'dynims'):
				print("Read in the dynamics first or to use T2w, set use_dyn to 0")
				return

			markupims=self.dynims[:,:,:,-1]

		if use_dyn==0:
			if not hasattr(self,'T2wims'):
				print("Read in the T2w data first - patient.read_T2w() or to use dynamics, set use_dyn to 1")
				return

			markupims=self.T2wims

		# Display the markup images one at a time and mark roi if required
		numslices=markupims.shape[2]
		mask=np.zeros((markupims.shape),dtype='uint16')

		if use_roipoly==1: #Method using roipoly for a detailed mask drawing
			for jj in range(numslices):
				h3=plt.figure()
				plt.imshow(markupims[:,:,jj],cmap='gray')
				useslice=input('Use this slice? 1=yes, 0=no')
				if int(useslice)==0:
					plt.close(h3)
				elif int(useslice)==1:
					roi=roipoly.roipoly()
					input('press enter to continue')
					mask[:,:,jj]=roi.getMask(markupims[:,:,jj])
				else:
					print('Aborting')
					return

		if use_roipoly==0: # Method using an ellipse to exclude unwanted bits of the image from calculation
			import roi_ellipse
			mask=roi_ellipse.draw_ellipse(markupims)



		if use_dyn==1 and tightmask==0:
			self.dynmask=mask
			maskname='dynmask.npy'
		if use_dyn==1 and tightmask==1:
			self.dyntightmask=mask
			maskname='dyntightmask'

		if use_dyn==0 and tightmask==0:
			self.mask=mask
			maskname='mask.npy'
		if use_dyn==0 and tightmask==1:
			self.tightmask=mask
			maskname='tightmask.npy'


		if save==1:
			print('Saving '+maskname)
			np.save(os.path.join(self.patientdirect,'Analysis',maskname),mask)


	def convert_mask(self,tightmask=0):
		#Method to convert the T2w mask to the size of the dynamic images
		if not hasattr(self,'dyninfo'):
			print('Need dynamic data to convert the mask')
			return
		
		if tightmask==1:
			if not hasattr(self,'tightmask'):
				print('need to make a tightmask first')
				return

		elif tightmask==0:
			if not hasattr(self,'mask'):
				print('No mask defined yet')
				return
		
		smallmask=np.zeros(self.dynims.shape[0:3],dtype='uint16')
		numslices=self.dynims.shape[2]
		if tightmask==0:
			largemask=self.mask
		if tightmask==1:
			largemask=self.tightmask

		for kk in range(numslices):
			smallmask[:,:,kk]=scipy.misc.imresize(largemask[:,:,kk],self.dynims.shape[0:2],interp='nearest')/255

		if tightmask==0:
			self.dynmask=smallmask
			np.save(os.path.join(self.patientdirect,'Analysis','dynmask.npy'),smallmask)
		elif tightmask==1:
			self.dyntightmask=smallmask
			np.save(os.path.join(self.patientdirect,'Analysis','dyntightmask.npy'),smallmask)

	def load_mask(self):
		#method to load previously made dynamic mask
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','mask.npy')):
			print('No mask saved')
			return
		else:
			self.mask=np.load(os.path.join(self.patientdirect,'Analysis','mask.npy'))

	def load_dynmask(self):
		#method to load previously made dynamic mask
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','dynmask.npy')):
			print('No mask saved')
			return
		else:
			self.dynmask=np.load(os.path.join(self.patientdirect,'Analysis','dynmask.npy'))

	def load_tightmask(self):
		#method to load previously made tight dynamic mask
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','dyntightmask.npy')):
			print('No dynamic tight mask saved')
			return
		else:
			self.dyntightmask=np.load(os.path.join(self.patientdirect,'Analysis','dyntightmask.npy'))


	def read_analyze_roi(self, debug=0, UMCU=0):
		
		if not hasattr(self,'dyninfo'):
			print('Read in the dynamic data to check for ROI placement')
			return
		# Method to read in analyze img hdr rois from Sjogren's data
		# Find the associated text file for checking later
		#txtfile=glob.glob(os.path.join(self.patientdirect,'Analysis','ROI','*.txt'))[0]
		
		import ResampleImage
		import nibabel
		# Look in pt_directory/analysis/roi to find MPRAGE and SEG files from GT
		imgfiles=glob.glob(os.path.join(self.patientdirect,'Analysis','ROI','*.nii'))
		print('reading...')
		print('roi file: '+imgfiles[1])
		print('MPRAGE file: '+imgfiles[0])
		# Read the files
		roi=nibabel.load(imgfiles[1])
		MPRAGE=nibabel.load(imgfiles[0])
		# extract pixel values
		roidata=np.squeeze(roi.get_data())
		MPRAGEdata=np.squeeze(MPRAGE.get_data())
		if debug==1:
			print('MPRAGE from nifti shape: ')
			print(MPRAGEdata.shape)

		if glob.glob(os.path.join(self.dicomdirect,'*_MPRAGEp*')):
			MPRAGEfolder=glob.glob(os.path.join(self.dicomdirect,'*_MPRAGEp*'))[0]
			Dynfolder=glob.glob(os.path.join(self.dicomdirect,'*minutes*'))[0]
		else:
			MPRAGEfolder=glob.glob(os.path.join(self.dicomdirect,'*3D T1'))[0]
			Dynfolder=glob.glob(os.path.join(self.dicomdirect,'*dynamic'))[0]


		#Reorient data from nifti format IMANOVA:(sag rotated anticlockwise 90 deg, flip lr, slice dirn flipped round) to match dicom (sag)
		# UMCU: (ax rotated 90deg, flip lr) to match dicom (ax)
		print('Reorienting nifti data')
		print(MPRAGEdata.dtype)
		
		#MPRAGEdatanew=np.moveaxis(MPRAGEdata,0,2)
		MPRAGEdatanew=np.fliplr(np.rot90(MPRAGEdata,-1)).astype('<f4')
		if UMCU==0:
			MPRAGEdatanew=np.flip(MPRAGEdatanew,2)

		if debug==1:
			print('MPRAGE nifti datatype:')
			print(MPRAGEdatanew.dtype)
			print('MPRAGE nifti shape:')
			print(MPRAGEdatanew.shape)
			plt.figure()
			plt.gca().set_title('MPRAGE from nifti, reoriented')
			plt.imshow(MPRAGEdatanew[:,:,100])	
		
		#Reorient roi data from nifti
		#roidatanew=np.moveaxis(roidata,0,2)
		roidatanew=np.fliplr(np.rot90(roidata,-1))
		if UMCU==0:
			roidatanew=np.flip(roidatanew,2)

		#print(roidatanew.dtype)

		print('Interpolating MPRAGE...')
		# Interpolate MPRAGE DICOM image data with dynamic grid
		MPRAGEasDyn=ResampleImage.resample_dicom(Dynfolder,MPRAGEfolder,MPRAGEdatanew,self.dynims.shape[2],MPRAGEdatanew.shape[2])
		#MPRAGEasDyn=np.load(os.path.join(self.patientdirect,'Analysis','MPRAGEasDyn.npy'))

		print('Interpolating ROI...')
		# Interpolate MPRAGE roi data with dynamic grid
		ROIasDyn=ResampleImage.resample_dicom(Dynfolder,MPRAGEfolder,roidatanew,self.dynims.shape[2],MPRAGEdatanew.shape[2])
		#ROIasDyn=np.load(os.path.join(self.patientdirect,'Analysis','ROIasDyn.npy'))

		# Reorient the data to appear as we would expect (L on right of image)
		#MPRAGEasDyn=np.rot90(MPRAGEasDyn,1)
		#ROIasDyn=np.rot90(ROIasDyn,1)

		#Quick visual check on one slice, also print text file
		# import overlay
		if debug==1:
			plt.figure()
			if UMCU==0:
				overlay.overlay(MPRAGEdatanew[183,:,:],roidatanew[183,:,:],0.2,6)
			else:
				overlay.overlay(MPRAGEdatanew[:,:,50],roidatanew[:,:,50],0.2,6)
			plt.gca().set_title('MPRAGE and rois from nifti')
			plt.figure()
			overlay.overlay(MPRAGEasDyn[:,:,6],ROIasDyn[:,:,6],0.2,6)
			plt.gca().set_title('MPRAGE and rois from nifti, both resampled')
			plt.figure()
			overlay.overlay(self.dynims[:,:,6,10],ROIasDyn[:,:,6],0.2,6)
			plt.gca().set_title('Dynamic from dicom, roi resampled from nifti')

		# with open(txtfile) as f:
		# 	for line in f:
		# 		print(line)

		#ROIinfo=np.genfromtxt(txtfile,delimiter="'",usecols=1,dtype=str,autostrip='True')
		#ROIinfo=np.append(ROIinfo[1:],'parotid')
		#ROIinfo=np.append(ROIinfo,'submandibular')
		#self.ROIinfo=ROIinfo

		self.MPRAGEasDyn=MPRAGEasDyn
		self.ROIasDyn=ROIasDyn

		np.save(os.path.join(self.patientdirect,'Analysis','MPRAGEasDyn.npy'),MPRAGEasDyn)
		np.save(os.path.join(self.patientdirect,'Analysis','ROIasDyn.npy'),ROIasDyn)
		#np.save(os.path.join(self.patientdirect,'Analysis','ROIinfo.npy'),ROIinfo)


	def make_T1map(self,seqtype='IR',measnum=0, save=1): # method to do T1 fitting, seqtype is either IR, SR, IRTrueFISP or VFA
		#if not hasattr(self, 'T1data'):
		#	print('Read the T1 data first - patient.read_T1data()')
		#	return
		if measnum==0:
			T1info=self.T1info
		if measnum==1:
			T1info=self.T1info2

		TIs=T1info['TIs'][0]
		TR=T1info['TR'][0]
		N=T1info['N'][0]
		flip=T1info['FlipAngle']
		sequenceparams=(flip,np.ceil(N/2),TR,N,4000)
		
		map=np.zeros(self.dynmask.shape)
		#self.T1map=np.zeros(self.dynmask.shape)

		if measnum==0:
			T1data=self.T1data
		if measnum==1:
			T1data=self.T1data2

		data=np.zeros((self.dynmask.shape+(1,)))
		data[:,:,:,0]=self.dynmask
		data=np.tile(data,(1,1,1,len(TIs)))
		data=T1data*data

		for sl in range(T1data.shape[2]): #for each slice
			print(sl)
			if np.sum(self.dynmask[:,:,sl])!=0: #if there are pixels in the slice
				for i in range(T1data.shape[0]):
					for j in range(T1data.shape[1]):
						curve=np.squeeze(data[i,j,sl,:])
						if np.sum(curve)!=0:
							curve=curve/np.amax(curve)
							if seqtype=='IR':
								fit=IRTurboFLASH.fittingfun(TIs,sequenceparams,curve)
							if seqtype=='IRTrueFISP':
								fit=IRTrueFISP.fittingfun(TIs,TR,curve)
							if seqtype=='VFA':
								fit=FLASH.fittingfun(flip,TR,curve/np.max(curve),startguess=[10,500])
								# Swap M0 and T1 around to match other fitting algorithms
								tmp=fit.x[0]
								fit.x[0]=fit.x[1]
								fit.x[1]=tmp
							if seqtype=='SR':
								fit=SRTurboFLASH.fittingfun(TIs,TR,flip,np.ceil(N/2),curve)
							map[i,j,sl]=fit.x[0]
							
		# Save the map

		if measnum==0:
			self.T1map=map
			mapname='T1map.npy'
		if measnum==1:
			self.T1map2=map
			mapname='T1map2.npy'

		if save==1:
			np.save(os.path.join(self.patientdirect,'Analysis',mapname),map)



	def load_T1map(self,measnum=0):
		#method to load previously made T1 map
		if measnum==0:
			filename='T1map.npy'
		elif measnum==1:
			filename='T1map2.npy'

		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis',filename)):
			print('No T1 map saved')
			return
		
		if measnum==0:
			self.T1map=np.load(os.path.join(self.patientdirect,'Analysis',filename))
		if measnum==1:
			self.T1map2=np.load(os.path.join(self.patientdirect,'Analysis',filename))
	
	def get_iAUC(self, save=1, baselinepts=10):
		#method to calculate iAUC for the pixels within the mask, set save to 1 to save to a file
		if not hasattr(self,'T1map'):
			print('Calculate the T1 map first')
			return
		if not hasattr(self,'dynims'):
		 	print('Load the dynamic images first')
		 	return

	 	#find TR and flip angle
		flip=self.dyninfo['FlipAngle']
		TR=self.dyninfo['TR']/1000
		tres=self.dyninfo['tres'][0]
		if tres==0:
			print('Time resolution is currently zero - correct this before continuing')
			return

		iAUC=np.zeros(self.dynmask.shape)

		for sl in range(self.dynmask.shape[2]): #for each slice
			print(sl)
			if np.sum(self.dynmask[:,:,sl])!=0: #if there are pixels in the slice
				for i in range(self.dynmask.shape[0]): #loop over rows and cols
					for j in range(self.dynmask.shape[1]):
							if self.dynmask[i,j,sl]==1:
								uptakeConc=FLASH.SI2Conc(np.squeeze(self.dynims[i,j,sl,:]),TR,flip,self.T1map[i,j,sl]/1000,baselinepts,None)
								if np.sum(np.isnan(uptakeConc))==0 and np.sum(uptakeConc)>0:
									iAUC[i,j,sl]=np.trapz(uptakeConc,dx=tres)

			#plt.figure()
			#plt.imshow(iAUC[:,:,sl])

		# Save the map if requested:
		if save==1:
			np.save(os.path.join(self.patientdirect,'Analysis','iAUC.npy'),iAUC)
		
		self.iAUC=iAUC


	def load_iAUC(self):
		#method to load previously made iAUC map
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','iAUC.npy')):
			print('No iAUC saved')
			return
		else:
			self.iAUC=np.load(os.path.join(self.patientdirect,'Analysis','iAUC.npy'))

	def get_EnhancingFraction(self, save=1, baselinepts=10, threshold=0):

		#Method for enhancing fraction calculation
		if not hasattr(self,'dynims'):
			print("Use read_dynamics to load the dynamic images first")
			return
		if not hasattr(self,'dyntightmask'):
			print("make a tight mask over the tumour first, and convert to dynamic size")
			return
		tres=self.dyninfo['tres'][0]
		if tres==0:
			print('Time resolution is currently zero - correct this before continuing')
			return


		#Find the sd in the baseline
		baselinemean=np.mean(self.dynims[:,:,:,1:baselinepts],3)
		baselinestd=np.std(self.dynims[:,:,:,1:baselinepts],3)
		#plt.figure()
		#plt.imshow(baselinestd[:,:,10]*self.dyntightmask[:,:,10],interpolation='nearest')
		#plt.colorbar()

		#Find the mean std over the mask pixels
		numtightmaskpixels=np.sum(self.dyntightmask)
		meanbaselinestd=np.sum(baselinestd*self.dyntightmask)/numtightmaskpixels
		print('pixels in mask = '+str(numtightmaskpixels))
		print('baseline sd = '+str(meanbaselinestd))
		
		#Calculate threshold
		if threshold==0:
			threshold=meanbaselinestd*25
		print('threshold = '+str(threshold))

		#Find the maximum enhancement then mask to tightmask
		maxenhancement=np.ndarray.max(self.dynims,3)
		maxenhancement=maxenhancement*(self.dyntightmask) #maximum enhancement in SI
		#plt.figure()
		#plt.imshow(maxenhancement[:,:,10],vmax=200,interpolation='nearest')
		#plt.colorbar()

		#Threshold on mean baseline *3
		numenhancingpixels=np.sum(maxenhancement>threshold)
		print('number of enhancing pixels = '+str(numenhancingpixels))
		EnhancingFraction=numenhancingpixels/numtightmaskpixels
		print('Enhancing Fraction = '+str(EnhancingFraction))

		# Calculate the max enhancement as a fraction of baseline for comparison with other datasets
		maxenhancement=(maxenhancement/(baselinemean*self.dyntightmask))-1 #max enhancement over baseline (0 at start)

		# Calculate the enhancement as a fraction of the baseline (so, 0 at start) during the dynamic series
		percentenh=(self.dynims/np.repeat(baselinemean[:,:,:,np.newaxis],self.dynims.shape[3],3))-1
		#Also calculate the initial rate of enhancement using sav_golay filter
		IRE=np.zeros(maxenhancement.shape)
		IRE=np.amax(scipy.signal.savgol_filter(percentenh,window_length=3, polyorder=1, deriv=1, axis=3, delta=tres),axis=3)

		if save==1:
			np.save(os.path.join(self.patientdirect,'Analysis','MaxEnhancement.npy'),maxenhancement)
			np.save(os.path.join(self.patientdirect,'Analysis','EnhancingFraction.npy'),EnhancingFraction)
			np.save(os.path.join(self.patientdirect,'Analysis','InitialRateEnhancement.npy'),IRE)
		
		self.EnhancingFraction=EnhancingFraction
		self.MaxEnhancement=maxenhancement
		self.InitialRateEnhancement=IRE

	def get_EnhancingFractionConc(self, save=1, baselinepts=10, threshold=0):
		#Method for enhancing fraction calculation but with concentration rather than SI
		if not hasattr(self,'dynims'):
			self.read_dynamics('')
			
		if not hasattr(self,'dyntightmask'):
			self.load_tightmask()
			
		if not hasattr(self,'T1map'):
			self.load_T1map()
			

		#Make an array to hold the concentration curves:

		Conccurves=np.zeros(self.dynims.shape)

		#find TR and flip angle
		flip=self.dyninfo['FlipAngle']
		TR=self.dyninfo['TR']/1000
		tres=self.dyninfo['tres'][0]
		if tres==0:
			print('Time resolution is currently zero - correct this before continuing')
			return

		
		#Loop through to convert to concentration
		for sl in range(self.dyntightmask.shape[2]): #for each slice
			#print(sl)
			if np.sum(self.dyntightmask[:,:,sl])!=0: #if there are pixels in the slice
				for i in range(self.dyntightmask.shape[0]): #loop over rows and cols
					for j in range(self.dyntightmask.shape[1]):
							if self.dyntightmask[i,j,sl]==1:
								uptakeConc=FLASH.SI2Conc(np.squeeze(self.dynims[i,j,sl,:]),TR,flip,self.T1map[i,j,sl]/1000,baselinepts,None)
								if np.sum(np.isnan(uptakeConc))==0 and np.sum(uptakeConc)>0:
									Conccurves[i,j,sl,:]=uptakeConc

		
		#Find the sd in the baseline
		baselinestd=np.std(Conccurves[:,:,:,1:baselinepts],3)

		#Find the mean std over the mask pixels
		numtightmaskpixels=np.sum(self.dyntightmask)
		meanbaselinestd=np.sum(baselinestd*self.dyntightmask)/numtightmaskpixels
		print('pixels in mask = '+str(numtightmaskpixels))
		print('baseline sd = '+str(meanbaselinestd))
		#Calculate threshold
		if threshold==0:
			threshold=meanbaselinestd*25
		print('threshold = '+str(threshold))

		#Find the maximum enhancement then mask to tightmask
		maxenhancement=np.ndarray.max(Conccurves,3)
		maxenhancement=maxenhancement*self.dyntightmask
		# plt.figure()
		# plt.imshow(maxenhancement[:,:,10],interpolation='nearest',vmax=10)
		# plt.colorbar()

		numenhancingpixels=np.sum(maxenhancement>threshold)
		print('number of enhancing pixels = '+str(numenhancingpixels))
		EnhancingFraction=numenhancingpixels/numtightmaskpixels
		print('Enhancing Fraction Conc= '+str(EnhancingFraction))

		#Also calculate the initial rate of enhancement using sav_golay filter
		IRE=np.zeros(maxenhancement.shape)
		IRE=np.nanmax(scipy.signal.savgol_filter(Conccurves,window_length=3, polyorder=1, deriv=1, axis=3,delta=tres),axis=3)
		
		self.EnhancingFractionConc=EnhancingFraction
		self.MaxEnhancementConc=maxenhancement
		self.InitialRateEnhancementConc=IRE
		self.Conccurves=Conccurves

		if save==1:
			np.save(os.path.join(self.patientdirect,'Analysis','MaxEnhancementConc.npy'),maxenhancement)
			np.save(os.path.join(self.patientdirect,'Analysis','EnhancingFractionConc.npy'),EnhancingFraction)
			np.save(os.path.join(self.patientdirect,'Analysis','InitialRateEnhancementConc.npy'),IRE)


	def get_region_stats(self,image,mask): #Use a mask and an image to get mean, sd, median, quartiles and volume
		#mask=getattr(self,mask)
		#image=getattr(self,image)
		# Check image and mask are the same shape
		if mask.shape != image.shape:
			print('Image and mask need to be the same shape and size')
			return
		values_in_mask=image[mask==1]
		# Num pixels
		numpix=np.sum(mask)
		numnans=np.sum(np.isnan(values_in_mask))
		numinf=np.sum(np.isinf(values_in_mask))
		if numinf>0:
			values_in_mask=np.ma.masked_invalid(values_in_mask)
		# Mean
		regionmean=np.nanmean(values_in_mask)
		# SD
		regionsd=np.nanstd(values_in_mask)
		# Median
		regionmedian=np.nanmedian(values_in_mask)
		# LQ
		regionlq=np.nanpercentile(values_in_mask,25)
		# UQ
		regionuq=np.nanpercentile(values_in_mask,75)

		#print([numpix,numnans,regionmean,regionsd,regionmedian,regionlq,regionuq])
		self.get_region_stats_result=[numpix,numnans,regionmean,regionsd,regionmedian,regionlq,regionuq]
		

	def plot_mean_curve(self, image, mask, x=0): #Use a mask and an image to plot a mean curve from a region mask. Put time in x variable if required
		image=getattr(self,image)
		# Check image and mask are the same shape
		if mask.shape[0:2] != image.shape[0:2]:
			print('Image and mask need to be the same size in first 2 dimensions')
			return
		numpix=np.sum(mask)
		curve=np.mean(image[mask==1,:],0)
		if x==0: # If x was zero, just plot curve
			plt.plot(curve)
		if x!=0: # If not, use the x provided
			plt.plot(x,curve)
		self.meancurve=curve
		np.save(os.path.join(self.patientdirect,'Analysis','meancurve.npy'),curve)


	# Fitting
	#####################################################
	def fit_Kety(self):
		pass

	def fit_ExtKety(self, SIflag=0, save=1):
		# To fit SI, set SIflag to 1, and to use second T1 measurement to correct SI curves, set use2ndT1=1

		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF available')
			return

		# Check for T1 map(s)
		if not hasattr(self,'T1map'):
			print('No T1map available')
			return

		if not hasattr(self,'dynmask'):
			print('No mask available')
			return

		self.ExtKetyfitConc=np.zeros((self.dynmask.shape+(4,)))
		self.ExtKetyfitSI=np.zeros((self.dynmask.shape+(4,)))
		self.ExtKetyfitConcEnhFlag=np.zeros((self.dynmask.shape))
		self.ExtKetyfitSIEnhFlag=np.zeros((self.dynmask.shape))

		TR=self.dyninfo['TR']/1000
		flip=self.dyninfo['FlipAngle']
		

		if SIflag==1:
			for sl in range(self.dynims.shape[2]):
					count=0
					total=np.sum(self.dynmask[:,:,sl])
					for i in range(self.dynims.shape[0]):
						for j in range(0,self.dynims.shape[1]):
							if self.dynmask[i,j,sl]==1:
								count=count+1
								print('Pixel '+str(count)+' of '+str(total)+' in slice '+str(sl))
								uptake=np.squeeze(self.dynims[i,j,sl,:])
								T1base=self.T1map[i,j,sl]
								fit, chi2=ExtKety.ExtKetyfittingSI(self.t, self.AIF/(1-self.hct), uptake, 15, TR, flip, T1base/1000)
								self.ExtKetyfitSIEnhFlag[i,j,sl]=self.enhancement_flag(uptake, 3, chi2)
								self.ExtKetyfitSI[i,j,sl,:]=fit
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','ExtKetyfitSImaps.npy'),self.ExtKetyfitSI)
				np.save(os.path.join(self.patientdirect,'Analysis','ExtKetyfitSIEnhFlag.npy'),self.ExtKetyfitSIEnhFlag)


		else:
			for sl in range(self.dynims.shape[2]):
					for i in range(self.dynims.shape[0]):
						for j in range(0,self.dynims.shape[1]):
							if self.dynmask[i,j,sl]==1:
								uptakeConc=FLASH.SI2Conc(self.dynims[i,j,sl,:],TR,flip,self.T1map[i,j,sl]/1000,15,None)
								print(i,j,sl)
								if np.isnan(np.sum(uptakeConc))==0:
									ExtKetyfitConc, chi2=ExtKety.ExtKetyfittingConc(self.t, self.AIF/(1-self.hct), uptakeConc)
									self.ExtKetyfitConcEnhFlag[i,j,sl]=self.enhancement_flag(uptakeConc, 3, chi2)
									self.ExtKetyfitConc[i,j,sl,:]=ExtKetyfitConc

			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','ExtKetyfitConcmaps.npy'),self.ExtKetyfitConc)
				np.save(os.path.join(self.patientdirect,'Analysis','ExtKetyfitConcEnhFlag.npy'),self.ExtKetyfitConcEnhFlag)


	def fit_2CXM(self,SIflag,save,use2ndT1=0):
		# To fit SI, set SIflag to 1, and to use second T1 measurement to correct SI curves, set use2ndT1=1
		if SIflag==0 and use2ndT1==1:
			print('Concentration fitting with second T1 not yet implemented')
			return

		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF available')
			return

		# Check for T1 map(s)
		if not hasattr(self,'T1map'):
			print('No T1map available')
			return
		if use2ndT1==1 and not hasattr(self,'T1map2'):
			print('Second T1map not found')
			return

		if not hasattr(self,'dynmask'):
			print('No mask available')
			return

		self.TwoCXMfitConc=np.zeros((self.dynmask.shape+(7,)))
		self.TwoCXMfitSI=np.zeros((self.dynmask.shape+(7,)))
		
		TR=self.dyninfo['TR']/1000
		flip=self.dyninfo['FlipAngle']
		rflip=flip*np.pi/180
	
		if SIflag==1:
			# Check if second T1 map is to be used for correction, if so calculate correction matrix
			if use2ndT1==1:
				SI0=np.mean(self.dynims[:,:,:,1:15],3) #find mean of beginning
				SI1=np.mean(self.dynims[:,:,:,-15:],3) #find mean of end
				M0=FLASH.CalcM0(SI0,TR,flip,self.T1map/1000)
				SI1C=FLASH.SIeqn([M0,self.T1map2/1000],flip,TR) #Find corrected SI at the end of the dynamic from the post contrast T1 measurement
				k=(SI1C-SI0)/(SI1-SI0)

			# Now work through slices, rows and columns
			for sl in range(self.dynims.shape[2]):
				count=0
				total=np.sum(self.dynmask[:,:,sl])
				for i in range(self.dynims.shape[0]):
					for j in range(0,self.dynims.shape[1]):
						if self.dynmask[i,j,sl]==1:
							count=count+1
							print('Pixel '+str(count)+' of '+str(total)+' in slice '+str(sl))
							uptake=np.squeeze(self.dynims[i,j,sl,:])
							#if correcting with second T1map, do that now
							if use2ndT1==1:
								uptake=k[i,j,sl]*(uptake-SI0[i,j,sl])+SI0[i,j,sl]
							T1base=self.T1map[i,j,sl]
							fit=TwoCXM.TwoCXMfittingSI(self.t, self.AIF/(1-self.hct), uptake, None, 15, TR, flip, T1base/1000,[0.03,0.3])
							self.TwoCXMfitSI[i,j,sl,:]=fit
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','TwoCXMfitSImaps.npy'),self.TwoCXMfitSI)

		else:
			for sl in range(self.dynims.shape[2]):
					for i in range(self.dynims.shape[0]):
						for j in range(0,self.dynims.shape[1]):
							if self.dynmask[i,j,sl]==1:
								uptakeConc=FLASH.SI2Conc(self.dynims[i,j,sl,:],TR,flip,self.T1map[i,j,sl]/1000,15,None)
								print(i,j,sl)
								if np.isnan(np.sum(uptakeConc))==0:
									TwoCXMfitConc=TwoCXM.TwoCXMfittingConc(self.t, self.AIF/(1-self.hct), uptakeConc, None)
									self.TwoCXMfitConc[i,j,sl,:]=TwoCXMfitConc
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','TwoCXMfitConcmaps.npy'),self.TwoCXMfitConc)

	def fit_2CUM(self,SIflag,save,usebothT1=0,use2ndT1=0):
		# To fit SI, set SIflag to 1, and to use second T1 measurement to correct SI curves, set use2ndT1=1
		if SIflag==0 and (use2ndT1==1 or usebothT1==1):
			print('Concentration fitting with second T1 not yet implemented')
			return
		# To fit SI, set SIflag to 1
		# To save maps, set save flag to 1
		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF available')
			return

		if not hasattr(self,'T1map'):
			print('No T1map available')
			return
		if (use2ndT1==1 or usebothT1==1) and not hasattr(self,'T1map2'):
			print('Second T1map not found')
			return

		if not hasattr(self,'dynmask'):
			print('No mask available')
			return

		self.TwoCUMfitConc=np.zeros((self.dynmask.shape+(6,)))
		self.TwoCUMfitSI=np.zeros((self.dynmask.shape+(6,)))
		
		TR=self.dyninfo['TR']/1000
		flip=self.dyninfo['FlipAngle']
		

		if SIflag==1:
			# Check if second T1 map is to be used for correction, if so calculate correction matrix
			if usebothT1==1:
				SI0=np.mean(self.dynims[:,:,:,1:15],3) #find mean of beginning
				SI1=np.mean(self.dynims[:,:,:,-15:],3) #find mean of end
				M0=FLASH.CalcM0(SI0,TR,flip,self.T1map/1000)
				SI1C=FLASH.SIeqn([M0,self.T1map2/1000],flip,TR) #Find corrected SI at the end of the dynamic from the post contrast T1 measurement
				k=(SI1C-SI0)/(SI1-SI0)

			for sl in range(self.dynims.shape[2]):
				count=0
				total=np.sum(self.dynmask[:,:,sl])
				for i in range(self.dynims.shape[0]):
					for j in range(0,self.dynims.shape[1]):
						if self.dynmask[i,j,sl]==1:
							count=count+1
							print('Pixel '+str(count)+' of '+str(total)+' in slice '+str(sl))
							uptake=np.squeeze(self.dynims[i,j,sl,:])
							#if correcting with second T1map, do that now
							if usebothT1==1:
								uptake=k[i,j,sl]*(uptake-SI0[i,j,sl])+SI0[i,j,sl]
							T1base=self.T1map[i,j,sl]
							if use2ndT1==1:
								T1base=self.T1map2[i,j,sl]
							fit=TwoCUM.TwoCUMfittingSI(self.t, self.AIF/(1-self.hct), uptake, None, 15, TR, flip, T1base/1000,[0.0003,0.03],use2ndT1)
							self.TwoCUMfitSI[i,j,sl,:]=fit
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','TwoCUMfitSImapsBothT1.npy'),self.TwoCUMfitSI)

		else:
			for sl in range(self.dynims.shape[2]):
					for i in range(self.dynims.shape[0]):
						for j in range(0,self.dynims.shape[1]):
							if self.dynmask[i,j,sl]==1:
								uptakeConc=FLASH.SI2Conc(self.dynims[i,j,sl,:],TR,flip,self.T1map[i,j,sl]/1000,15,None)
								print(i,j,sl)
								if np.isnan(np.sum(uptakeConc))==0:
									TwoCUMfitConc=TwoCUM.TwoCUMfittingConc(self.t, self.AIF/(1-self.hct), uptakeConc, None)
									self.TwoCUMfitConc[i,j,sl,:]=TwoCUMfitConc
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','TwoCUMfitConcmaps.npy'),self.TwoCUMfitConc)

	def fit_AATH(self,SIflag,save):
		# To fit SI, set SIflag to 1
		# To save the maps as a numpy array, set save to 1 
		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF - if reading from file, patient.read_AIF_fromfittingfile')
			return
		self.AATHfitConc=np.zeros((self.dynmask.shape+(6,)))
		self.AATHfitSI=np.zeros((self.dynmask.shape+(6,)))

		TR=self.dyninfo['TR']/1000
		flip=self.dyninfo['FlipAngle']
		
		if SIflag==1:
			for sl in range(self.dynims.shape[2]):
					count=0
					total=np.sum(self.dynmask[:,:,sl])
					for i in range(self.dynims.shape[0]):
						for j in range(0,self.dynims.shape[1]):
							if self.dynmask[i,j,sl]==1:
								count=count+1
								print('Pixel '+str(count)+' of '+str(total)+' in slice '+str(sl))
								uptake=np.squeeze(self.dynims[i,j,sl,:])
								T1base=self.T1map[i,j,sl]
								fit=AATH.AATHfittingSI(self.t, self.AIF/(1-self.hct), uptake, None, 15, TR, flip, T1base/1000)
								self.AATHfitSI[i,j,sl,:]=fit
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','AATHfitSImaps.npy'),self.AATHfitSI)


		else:
			for sl in range(self.dynims.shape[2]):
					for i in range(self.dynims.shape[0]):
						for j in range(0,self.dynims.shape[1]):
							if self.dynmask[i,j,sl]==1:
								uptakeConc=FLASH.SI2Conc(self.dynims[i,j,sl,:],TR,flip,self.T1map[i,j,sl]/1000,15,None)
								print(i,j,sl)
								if np.isnan(np.sum(uptakeConc))==0:
									AATHfitConc=AATH.AATHfittingConc(self.t, self.AIF/(1-self.hct), uptakeConc, None)
									self.AATHfitConc[i,j,sl,:]=AATHfitConc

			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','AATHfitConcmaps.npy'),self.AATHfitConc)

	def enhancement_flag(self, curve, k, chi2_model):
		# Helper method to fit a constant to a curve and output whether this is a better fit, according to the Akaike criterion
		#Fit constant
		fit_const=np.polyfit(self.t,curve,0,full=True)
		chi2_const=fit_const[1][0]
		#print(fit_const, chi2_const)
		#print(fit_const[0][0])
		
		# Akaike = 2*k + n*ln(chi2)
		Akaike_const=(2) + len(curve)*np.log(chi2_const) #k=1 for a constant
		#print(Akaike_const)
		Akaike_model=(2*k) + len(curve)*np.log(chi2_model**2)
		#print(Akaike_model)
		# Return 1 if model is a better fit than constant, i.e if enhancing, 0 if constant is a better fit (i.e. not enhancing)
		flag=1
		#print(Akaike_const<Akaike_model)
		#print(fit_const[0][0]<0)
		if (Akaike_const<Akaike_model) or (fit_const[0][0]<0):
			flag=0

		return flag 


	def load_AATHConcparammaps(self):
		#Method to load maps from numpy arrays
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','AATHfitConcmaps.npy')):
			print('No AATH maps saved')
			return
		else:
			self.AATHfitConc=np.load(os.path.join(self.patientdirect,'Analysis','AATHfitConcmaps.npy'))

	def load_AATHSIparammaps(self):
		#Method to load maps from numpy arrays
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','AATHfitSImaps.npy')):
			print('No AATH maps saved')
			return
		else:
			self.AATHfitSI=np.load(os.path.join(self.patientdirect,'Analysis','AATHfitSImaps.npy'))

	def load_TwoCXMSIparammaps(self):
		#Method to load maps from numpy arrays
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','TwoCXMfitSImaps.npy')):
			print('No TwoCXM maps saved')
			return
		else:
			self.TwoCXMfitSI=np.load(os.path.join(self.patientdirect,'Analysis','TwoCXMfitSImaps.npy'))

	def load_TwoCUMSIparammaps(self):
		#Method to load maps from numpy arrays
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','TwoCUMfitSImaps.npy')):
			print('No TwoCUMSI maps saved')
			return
		else:
			self.TwoCUMfitSI=np.load(os.path.join(self.patientdirect,'Analysis','TwoCUMfitSImaps.npy'))

	def load_TwoCUMSIparammapsBothT1(self):
		#Method to load maps from numpy arrays
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','TwoCUMfitSImapsBothT1.npy')):
			print('No TwoCUMSIBothT1 maps saved')
			return
		else:
			self.TwoCUMfitSIBothT1=np.load(os.path.join(self.patientdirect,'Analysis','TwoCUMfitSImapsBothT1.npy'))

	def load_TwoCUMConcparammaps(self):
		#Method to load maps from numpy arrays
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','TwoCUMfitConcmaps.npy')):
			print('No TwoCUMConc maps saved')
			return
		else:
			self.TwoCUMfitConc=np.load(os.path.join(self.patientdirect,'Analysis','TwoCUMfitConcmaps.npy'))

	def write_iAUC(self, dynfoldertag):
		import time
		# Method to write iAUC to a dicom files, setting tags appropriately
		# Basetags folder tag is the tag for the dynamic series from which we will borrow slice positions etc
		
		DynFolder=glob.glob(os.path.join(self.dicomdirect,dynfoldertag))
		print('Borrowing dicom headers from here:')
		print(os.path.join(self.dicomdirect,DynFolder[0]))
		os.chdir(DynFolder[0])
		
		#Borrowing headers from the first dynamic volume
		filenames=glob.glob('*.dcm') # Find all dicom files in the folder
		iAUC=self.iAUC # Get the iAUC map to be written
		numslices=iAUC.shape[2] # Find the number of slices in the map

		#Read first file and check for Siemens or Philips
		tmp=dicom.read_file(filenames[0])

		if tmp[0x08,0x70].value=='SIEMENS': # If SIEMENS, use the first (numslices) files, this will be the first dynamic volume
			filenames=filenames[0:numslices]
		elif tmp[0x08,0x70].value=='Philips Medical Systems': # If Philips, find time points and use this to extract one dynamic volume
			timepoints=int(len(filenames)/numslices)
			filenames=filenames[0::timepoints]
		else:
			print('Unknown manufacturer')
			return

		#Make new series UID, creation time and date)
		Seriestime=time.strftime('%H%M%S')
		Seriesdate=time.strftime('%Y%m%d')
		SeriesUID=dicom.UID.generate_uid()

		for i in range(numslices):
			#Read original dicom
			tmp=dicom.read_file(filenames[i])
			
			#Prepare iAUC map for this slice
			map=iAUC[:,:,i]
			map=abs(map)
			#Remove values above 99th centile
			mapupperlim=np.percentile(map[map>0],99)
			map[map>mapupperlim]=mapupperlim
			map=map.astype('uint16')
			tmp.PixelData=map.tostring()

			# Set new header values
			tmp[0x08,0x08].value='DERIVED' #Image Type
			tmp[0x08,0x12].value=Seriesdate #Creation date
			tmp[0x08,0x13].value=Seriestime #Creation time
			tmp[0x08,0x18].value=dicom.UID.generate_uid() #Instance UID
			tmp[0x08,0x70].value='' # Manufacturer name
			tmp[0x08,0x1090].value='' # Manufacturer model
			tmp[0x08,0x103E].value='iAUC_map' # Description
			tmp[0x20,0x0e].value=SeriesUID # Series UID
			tmp[0x20,0x11].value=str(tmp[0x20,0x11].value)+'001' # Series number
			if tmp[0x08,0x70]=='SIEMENS':
				tmp[0x28,0x0106].value=np.min(map) # Maximum pixel value, set for Siemens only
				tmp[0x28,0x0107].value=np.max(map) # Minimum pixel value, set for Siemens only


			#Write the file
			dicom.write_file(str(tmp[0x20,0x11].value)+'_iAUC_'+str(i+1)+'.dcm',tmp)

		# Make a new folder and put these files into it
		newfiles=glob.glob('*iAUC*')
		newfoldername=str(tmp[0x20,0x11].value)+'_iAUC'
		print(newfoldername)
		os.mkdir(newfoldername)
		for x in newfiles:
			os.rename(x,os.path.join(newfoldername,x))
		shutil.move(newfoldername,'..')

	def write_maps(self,dynfoldertag,mapname,seriesnumstart=1):
		import time
		# Method to save specified maps (i.e. AATH, 2CXM, conc, SI etc) to dicom
		DynFolder=glob.glob(os.path.join(self.dicomdirect,dynfoldertag))		
		print('Borrowing dicom headers from here:')
		print(os.path.join(self.dicomdirect,DynFolder[0]))
		os.chdir(DynFolder[0])
		
		#Borrowing headers from the first dynamic volume
		filenames=glob.glob('*.dcm') # Find all dicom files in the folder
		maps=getattr(self,mapname)
		numslices=maps.shape[2] # Find the number of slices in the map
		nummaps=maps.shape[-1] # Find the number of maps to write out
		#7 maps: E, Fp, ve, vp, chi2, toff, status
		#6 maps: E, Fp, vp, chi2, toff, status

		#Read first file and check for Siemens or Philips
		tmp=dicom.read_file(filenames[0])

		if tmp[0x08,0x70].value=='SIEMENS': # If SIEMENS, use the first (numslices) files, this will be the first dynamic volume
			filenames=filenames[0:numslices]
		elif tmp[0x08,0x70].value=='Philips Medical Systems': # If Philips, find time points and use this to extract one dynamic volume
			timepoints=int(len(filenames)/numslices)
			filenames=filenames[0::timepoints]
		else:
			print('Unknown manufacturer')
			return

		#Get original series number and new times and Series UIDs for each map
		SeriesUIDs=[None]*6
		SeriesUIDs[0]=dicom.UID.generate_uid() #Must be done in separate calls with a pause, or it's all the same one!
		time.sleep(3)
		SeriesUIDs[1]=dicom.UID.generate_uid()
		time.sleep(2)
		SeriesUIDs[2]=dicom.UID.generate_uid()
		time.sleep(3)
		SeriesUIDs[3]=dicom.UID.generate_uid()
		time.sleep(2)
		SeriesUIDs[4]=dicom.UID.generate_uid()
		time.sleep(3)
		SeriesUIDs[5]=dicom.UID.generate_uid()
		#print(SeriesUIDs)

		Seriesdate=time.strftime('%Y%m%d')
		Seriestime=time.strftime('%H%M%S')
		OriginalSeriesNum=str(tmp[0x20,0x11].value)
		#print(OriginalSeriesNum)

		for i in range(numslices):
			#Read original dicom and set new header values that apply to all maps
			tmp=dicom.read_file(filenames[i])
			tmp[0x08,0x08].value='DERIVED' #Image Type		
			tmp[0x08,0x70].value='' # Manufacturer name
			tmp[0x08,0x1090].value='' # Manufacturer model

			#Prepare maps for this slice
			#E
			Emap=abs(maps[:,:,i,0])*100 #E in %
			Emap[Emap>150]=150 # Remove anything above 150
			Emap=Emap.astype('uint16')
			#Set the pixel data and remaining header values
			tmp.PixelData=Emap.tostring()
			tmp[0x08,0x12].value=Seriesdate #Creation date
			tmp[0x08,0x13].value=Seriestime #Creation time	
			tmp[0x08,0x18].value=dicom.UID.generate_uid() #Instance UID
			tmp[0x08,0x103E].value=mapname+'E_map' # Description
			tmp[0x20,0x0e].value=SeriesUIDs[0] # Series UID
			tmp[0x20,0x11].value=OriginalSeriesNum+np.str(seriesnumstart) # Series number
			if tmp[0x08,0x70].value=='SIEMENS':
				tmp[0x28,0x0106].value=np.min(map) # Maximum pixel value, set for Siemens only
				tmp[0x28,0x0107].value=np.max(map) # Minimum pixel value, set for Siemens only
			#Write the file
			dicom.write_file(str(tmp[0x20,0x11].value)+mapname+'_E_'+str(i+1)+'.dcm',tmp)

			#Fp
			Fpmap=abs(maps[:,:,i,1])*100*60*100 #Fp in ml/100ml/min multiplied by 100
			if np.sum(Fpmap>0):
				Fpmapupperlim=np.percentile(Fpmap[Fpmap>0],99) #Remove values above 99th centile
				Fpmap[Fpmap>Fpmapupperlim]=Fpmapupperlim
			Fpmap=Fpmap.astype('uint16')
			#Set the pixel data and remaining header values
			tmp.PixelData=Fpmap.tostring()
			tmp[0x08,0x12].value=Seriesdate #Creation date
			tmp[0x08,0x13].value=Seriestime #Creation time	
			tmp[0x08,0x18].value=dicom.UID.generate_uid() #Instance UID
			tmp[0x08,0x103E].value=mapname+'Fp_map' # Description
			tmp[0x20,0x0e].value=SeriesUIDs[1] # Series UID
			tmp[0x20,0x11].value=OriginalSeriesNum+np.str(seriesnumstart+1) # Series number
			if tmp[0x08,0x70].value=='SIEMENS':
				tmp[0x28,0x0106].value=np.min(map) # Maximum pixel value, set for Siemens only
				tmp[0x28,0x0107].value=np.max(map) # Minimum pixel value, set for Siemens only
			#Write the file
			dicom.write_file(str(tmp[0x20,0x11].value)+mapname+'_Fp_'+str(i+1)+'.dcm',tmp)

			#ve
			if nummaps==7:
				vemap=abs(maps[:,:,i,2])*100 #ve in %
				vemap[vemap>150]=150 #Remove anything above 150
				vemap=vemap.astype('uint16')
				#Set the pixel data and remaining header values
				tmp.PixelData=vemap.tostring()
				tmp[0x08,0x12].value=Seriesdate #Creation date
				tmp[0x08,0x13].value=Seriestime #Creation time	
				tmp[0x08,0x18].value=dicom.UID.generate_uid() #Instance UID
				tmp[0x08,0x103E].value=mapname+'ve_map' # Description
				tmp[0x20,0x0e].value=SeriesUIDs[2] # Series UID
				tmp[0x20,0x11].value=OriginalSeriesNum+np.str(seriesnumstart+2) # Series number
				if tmp[0x08,0x70].value=='SIEMENS':
					tmp[0x28,0x0106].value=np.min(map) # Maximum pixel value, set for Siemens only
					tmp[0x28,0x0107].value=np.max(map) # Minimum pixel value, set for Siemens only		
				#Write the file
				dicom.write_file(str(tmp[0x20,0x11].value)+mapname+'_ve_'+str(i+1)+'.dcm',tmp)

			#vp
			vpmap=abs(maps[:,:,i,-4])*100 #vp in %
			vpmap[vpmap>150]=150 #Remove anything above 150
			vpmap=vpmap.astype('uint16')
			#Set the pixel data and remaining header values
			tmp.PixelData=vpmap.tostring()
			tmp[0x08,0x12].value=Seriesdate #Creation date
			tmp[0x08,0x13].value=Seriestime #Creation time	
			tmp[0x08,0x18].value=dicom.UID.generate_uid() #Instance UID
			tmp[0x08,0x103E].value=mapname+'vp_map' # Description
			tmp[0x20,0x0e].value=SeriesUIDs[-3] # Series UID
			tmp[0x20,0x11].value=OriginalSeriesNum+np.str(seriesnumstart+3) # Series number
			if tmp[0x08,0x70].value=='SIEMENS':
				tmp[0x28,0x0106].value=np.min(map) # Maximum pixel value, set for Siemens only
				tmp[0x28,0x0107].value=np.max(map) # Minimum pixel value, set for Siemens only		
			#Write the file
			dicom.write_file(str(tmp[0x20,0x11].value)+mapname+'_vp_'+str(i+1)+'.dcm',tmp)	

			#toff
			toffmap=abs(maps[:,:,i,-2])*1000 #toff in milliseconds (max was 20 or 30 s)
			toffmap[toffmap>30000]=30000 #Remove anything above 30000 ms
			toffmap=toffmap.astype('uint16')
			#Set the pixel data and remaining header values
			tmp.PixelData=toffmap.tostring()
			tmp[0x08,0x12].value=Seriesdate #Creation date
			tmp[0x08,0x13].value=Seriestime #Creation time	
			tmp[0x08,0x18].value=dicom.UID.generate_uid() #Instance UID
			tmp[0x08,0x103E].value=mapname+'toff_map' # Description
			tmp[0x20,0x0e].value=SeriesUIDs[-1] # Series UID
			tmp[0x20,0x11].value=OriginalSeriesNum+np.str(seriesnumstart+4) # Series number
			if tmp[0x08,0x70].value=='SIEMENS':
				tmp[0x28,0x0106].value=np.min(map) # Maximum pixel value, set for Siemens only
				tmp[0x28,0x0107].value=np.max(map) # Minimum pixel value, set for Siemens only				
			#Write the file
			dicom.write_file(str(tmp[0x20,0x11].value)+mapname+'_toff_'+str(i+1)+'.dcm',tmp)

			#PS
			PSmap=abs(-1*np.log(1-maps[:,:,i,0])*maps[:,:,i,1])*60*100*100 #PS in 100*ml/100ml/min
			if np.sum(PSmap>0):
				PSmapupperlim=np.percentile(PSmap[PSmap>0],99) #Remove values above 99th percentile
				PSmap[PSmap>PSmapupperlim]=PSmapupperlim
			PSmap=PSmap.astype('uint16')
			#Set the pixel data and remaining header values
			tmp.PixelData=PSmap.tostring()
			tmp[0x08,0x12].value=Seriesdate #Creation date
			tmp[0x08,0x13].value=Seriestime #Creation time	
			tmp[0x08,0x18].value=dicom.UID.generate_uid() #Instance UID
			tmp[0x08,0x103E].value=mapname+'PS_map' # Description
			tmp[0x20,0x0e].value=SeriesUIDs[-2] # Series UID
			tmp[0x20,0x11].value=OriginalSeriesNum+np.str(seriesnumstart+5) # Series number
			if tmp[0x08,0x70].value=='SIEMENS':
				tmp[0x28,0x0106].value=np.min(map) # Maximum pixel value, set for Siemens only
				tmp[0x28,0x0107].value=np.max(map) # Minimum pixel value, set for Siemens only	
			#Write the file
			dicom.write_file(str(tmp[0x20,0x11].value)+mapname+'_PS_'+str(i+1)+'.dcm',tmp)

		dicomdirect=self.patientdirect+'/DICOM'
		newfiles=glob.glob('*_E_*')
		dirname=mapname+'_E'
		os.mkdir(dirname)
		for x in newfiles:
			os.rename(x,os.path.join(dirname,x))	
		
		print(os.path.exists(dicomdirect+'/'+dirname))	
		if os.path.exists(dicomdirect+'/'+dirname):
			print('Removing old maps')
			shutil.rmtree(os.path.join(dicomdirect,dirname))
		shutil.move(dirname,'..')

		newfiles=glob.glob('*_Fp_*')
		dirname=mapname+'_Fp'
		os.mkdir(dirname)
		for x in newfiles:
			os.rename(x,os.path.join(dirname,x))
		if os.path.exists(dicomdirect+'/'+dirname):
			shutil.rmtree(os.path.join(dicomdirect,dirname))
		shutil.move(dirname,'..')

		if nummaps==7:
			newfiles=glob.glob('*_ve_*')
			dirname=mapname+'_ve'
			os.mkdir(dirname)
			for x in newfiles:
				os.rename(x,os.path.join(dirname,x))
			if os.path.exists(dicomdirect+'/'+dirname):
				shutil.rmtree(os.path.join(dicomdirect,dirname))
			shutil.move(dirname,'..')

		newfiles=glob.glob('*_vp_*')
		dirname=mapname+'_vp'
		os.mkdir(dirname)
		for x in newfiles:
			os.rename(x,os.path.join(dirname,x))
		if os.path.exists(dicomdirect+'/'+dirname):
			shutil.rmtree(os.path.join(dicomdirect,dirname))
		shutil.move(dirname,'..')
		
		newfiles=glob.glob('*_toff_*')
		dirname=mapname+'_toff'
		os.mkdir(dirname)
		for x in newfiles:
			os.rename(x,os.path.join(dirname,x))
		if os.path.exists(dicomdirect+'/'+dirname):
			shutil.rmtree(os.path.join(dicomdirect,dirname))
		shutil.move(dirname,'..')

		newfiles=glob.glob('*_PS_*')
		dirname=mapname+'_PS'
		os.mkdir(dirname)
		for x in newfiles:
			os.rename(x,os.path.join(dirname,x))
		if os.path.exists(dicomdirect+'/'+dirname):
			shutil.rmtree(os.path.join(dicomdirect,dirname))
		shutil.move(dirname,'..')
		


	def apply_analyze_rois(self, UMCU=0):	
		import Sjogrens_export
		#Method to apply the read in rois to the following maps:
		#T1
		#MaxEnh
		#MaxEnhConc
		#Initial rate of enhancement
		#Initial rate of enhancement Conc
		#Ktrans
		
		#Check if maps are read in, if not then read in if present

		if (not hasattr(self,'T1map')): # If the map isn't present....
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','T1map.npy'))): #...and the file exists...
				print('loading T1 map')
				self.load_T1map() #... load it
			# else: # Or report that it's not found
			# 	print('T1 map not found')

		if (not hasattr(self,'MaxEnhancement')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','MaxEnhancement.npy'))):
				print('Loading MaxEnh')
				self.MaxEnhancement=np.load(os.path.join(self.patientdirect,'Analysis','MaxEnhancement.npy'))
			# else:
			# 	print('Maximum Enhancement map not found')

		if (not hasattr(self,'MaxEnhancementConc')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','MaxEnhancementConc.npy'))):
				print('Loading MaxEnhConc')
				self.MaxEnhancementConc=np.load(os.path.join(self.patientdirect,'Analysis','MaxEnhancementConc.npy'))
			# else:
			# 	print('Maximum Enhancement Conc map not found')

		if (not hasattr(self,'InitialRateEnhancement')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','InitialRateEnhancement.npy'))):
				print('Loading InitialRate Enhancement')
				self.InitialRateEnhancement=np.load(os.path.join(self.patientdirect,'Analysis','InitialRateEnhancement.npy'))
			# else:
			# 	print('Initial Rate Enhancement map not found')

		if (not hasattr(self,'InitialRateEnhancementConc')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','InitialRateEnhancementConc.npy'))):
				print('Loading Initial Rate enhancement Conc')
				self.InitialRateEnhancementConc=np.load(os.path.join(self.patientdirect,'Analysis','InitialRateEnhancementConc.npy'))
			# else:
			# 	print('Initial Rate Enhancement Conc map not found')

		if (not hasattr(self,'ExtKetyfitSI')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','ExtKetyfitSImaps.npy'))):
				print('Loading Ext KetyfitSImaps')
				self.ExtKetyfitSI=np.load(os.path.join(self.patientdirect,'Analysis','ExtKetyfitSImaps.npy'))
			# else:
			# 	print('Ext Kety fit SI maps not found')

		if (not hasattr(self,'ExtKetyfitSIEnhFlag')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','ExtKetyfitSIEnhFlag.npy'))):
				print('Loading ExtKetyfitSIEnhFlag')
				self.ExtKetyfitSIEnhFlag=np.load(os.path.join(self.patientdirect,'Analysis','ExtKetyfitSIEnhFlag.npy'))
			# else:
			# 	print('Ext Kety fit SI maps not found')

		if (not hasattr(self,'ExtKetyfitConc')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','ExtKetyfitConcmaps.npy'))):
				print('Loading ExtKetyfitConc')
				self.ExtKetyfitConc=np.load(os.path.join(self.patientdirect,'Analysis','ExtKetyfitConcmaps.npy'))
			# else:
			# 	print('Ext Kety fit Conc maps not found')

		if (not hasattr(self,'ExtKetyfitConcEnhFlag')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','ExtKetyfitConcEnhFlag.npy'))):
				print('Loading ExtKetyfitConcEnhFlag')
				self.ExtKetyfitConcEnhFlag=np.load(os.path.join(self.patientdirect,'Analysis','ExtKetyfitConcEnhFlag.npy'))
			# else:
			# 	print('Ext Kety fit Conc maps not found')

		if (not hasattr(self,'ROIasDyn')):
			if (os.path.isfile(os.path.join(self.patientdirect,'Analysis','ROIasDyn.npy'))):
				print('Loading ROI file')
				self.ROIasDyn=np.load(os.path.join(self.patientdirect,'Analysis','ROIasDyn.npy'))
			# else:
			# 	print('Ext Kety fit Conc maps not found')

		self.ROIinfo=np.array(['rt_parotid', 'lt_parotid', 'rt_submandibular', 'lt_submandibular', 'parotid', 'submandibular'], dtype='<U16')

		SEQ=1 # SEQ needs to be incremented for each result written to file
		
		if UMCU==1:
			subjectID=self.patientdirect[split('t')][-1] # Extract subject ID to write to file
		else:
			subjectID=self.patientdirect.split('_')[-2][-3:] # Extract subject ID to write to file

		#Work through masks and maps, printing out results and also sending to a csv file formatted for the spreadsheet of doom
		for X in range(6): # six possible masks - 1, 2, 3, 4, 1&2, 3&4
			#Set mask according to which region we're going for
			if X==0:
				mask=scipy.ndimage.morphology.binary_erosion(self.ROIasDyn==1,structure=np.ones((3,3,1)))
				reducedmask=mask*self.ExtKetyfitConcEnhFlag
			if X==1:
				mask=scipy.ndimage.morphology.binary_erosion(self.ROIasDyn==2,structure=np.ones((3,3,1)))
				reducedmask=mask*self.ExtKetyfitConcEnhFlag
			if X==2:
				mask=scipy.ndimage.morphology.binary_erosion(self.ROIasDyn==3,structure=np.ones((3,3,1)))
				reducedmask=mask*self.ExtKetyfitConcEnhFlag
			if X==3:
				mask=scipy.ndimage.morphology.binary_erosion(self.ROIasDyn==4,structure=np.ones((3,3,1)))
				reducedmask=mask*self.ExtKetyfitConcEnhFlag
			if X==4:
				mask=scipy.ndimage.morphology.binary_erosion((self.ROIasDyn==1)+(self.ROIasDyn==2),structure=np.ones((3,3,1)))
				reducedmask=mask*self.ExtKetyfitConcEnhFlag
			if X==5:
				mask=scipy.ndimage.morphology.binary_erosion((self.ROIasDyn==3)+(self.ROIasDyn==4),structure=np.ones((3,3,1)))
				reducedmask=mask*self.ExtKetyfitConcEnhFlag


			if np.sum(mask)>0:
				if hasattr(self, 'T1map'):
					self.get_region_stats(self.T1map,mask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', T1map, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])
					Sjogrens_export.send_to_file(subjectID,np.array(self.get_region_stats_result)[2:],'T1map',str(self.ROIinfo[X]),SEQ)
					SEQ=SEQ+5
					self.get_region_stats(self.T1map,reducedmask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', T1map reduced, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])

				if hasattr(self, 'MaxEnhancement'):
					self.get_region_stats(self.MaxEnhancement,reducedmask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', MaxEnhancement, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])
					Sjogrens_export.send_to_file(subjectID,np.array(self.get_region_stats_result)[2:]*100,'MaxEnhancement',str(self.ROIinfo[X]),SEQ) #fractional multiplied by 100 to get %
					SEQ=SEQ+5
					self.get_region_stats(self.MaxEnhancement,mask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', MaxEnhancement (not reduced), '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])

				# if hasattr(self, 'MaxEnhancementConc'):
				# 	erroneousindices=self.MaxEnhancementConc==0
				# 	self.MaxEnhancementConc[erroneousindices]=np.nan
				# 	self.get_region_stats(self.MaxEnhancementConc,mask)
				# 	print(self.patientdirect.split('/')[-1]+', MaxEnhancementConc, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result))
				# 	Sjogrens_export.send_to_file(subjectID,np.array(self.get_region_stats_result)[2:],'MaxEnhancementConc',str(self.ROIinfo[X]),SEQ)# s^-1 / 3.6MMOL^-1 s^-1 = MMOL
				# 	SEQ=SEQ+5

				if hasattr(self, 'InitialRateEnhancement'):
					self.get_region_stats(self.InitialRateEnhancement,reducedmask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', InitialRateEnhancement, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])
					Sjogrens_export.send_to_file(subjectID,np.array(self.get_region_stats_result)[2:]*100,'InitialRateEnhancement',str(self.ROIinfo[X]),SEQ)#Fractional multiplied by 100 to get %/s
					SEQ=SEQ+5
					self.get_region_stats(self.InitialRateEnhancement, mask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', InitialRateEnhancement (not reduced), '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])

				# if hasattr(self, 'InitialRateEnhancementConc'):
				# 	erroneousindices=self.InitialRateEnhancementConc==0
				# 	self.InitialRateEnhancementConc[erroneousindices]=np.nan
				# 	self.get_region_stats(self.InitialRateEnhancementConc,mask)
				# 	print(self.patientdirect.split('/')[-1]+', InitialRateEnhancementConc, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result))
				# 	Sjogrens_export.send_to_file(subjectID,np.array(self.get_region_stats_result)[2:]/3.6,'InitialRateEnhancementConc',str(self.ROIinfo[X]),SEQ) # s^-2 / 3.6MMOL^-1 s^-1 = MMOL/s
				# 	SEQ=SEQ+5

				# if hasattr(self, 'ExtKeyfitSI'):
				# 	self.get_region_stats(self.ExtKetyfitSI[:,:,:,0],mask)
				# 	print(self.patientdirect.split('/')[-1]+', ExtKeyfitSI, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result))			
				#	Sjogrens_export.send_to_file(subjectID,np.array(self.get_region_stats_result)[2:]*60,'ExtKetyfitSI',str(self.ROIinfo[X]),SEQ)
				#	SEQ=SEQ+5

				if hasattr(self, 'ExtKetyfitConc'):
					erroneousindices=self.ExtKetyfitConc[:,:,:,0]<0.0000101
					self.ExtKetyfitConc[erroneousindices]=np.nan
					self.get_region_stats(self.ExtKetyfitConc[:,:,:,0],reducedmask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', ExtKetyfitConc, '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])
					Sjogrens_export.send_to_file(subjectID,np.array(self.get_region_stats_result)[2:]*60,'ExtKetyfitConc',str(self.ROIinfo[X]),SEQ) # s^-1 *60 = min^-1
					SEQ=SEQ+5
					self.get_region_stats(self.ExtKetyfitConc[:,:,:,0],mask)
					with open ('/Users/lkershaw/Desktop/newresultsextended.csv','a',newline='') as csvfileextended:
						extendedresultswriter=csv.writer(csvfileextended,delimiter=',')
						extendedresultswriter.writerow([self.patientdirect.split('/')[-1]+', ExtKetyfitConc (not reduced), '+str(self.ROIinfo[X])+', '+str(self.get_region_stats_result)])

	

	def fit_TH(self):
		pass







