# Main file for defining the patient class to do DCE processing pixel by pixel

# Imports
import numpy as np
import dicom
import matplotlib.pyplot as plt
import os
import glob
import shutil
from tkinter import Tk
from tkinter import filedialog
from roipoly import roipoly

import scipy.misc
import SRTurboFLASH
import IRTurboFLASH
import scipy.optimize
import TwoCXM
import TwoCUM
import AATH
import FLASH
import getAIF
import pprint

class patient(object): # patient inherits from the object class
	
	def __init__(self,patientdirectory='select'):
		self.patientdirect=patientdirectory
		self.dicomdirect=os.path.join(self.patientdirect,'DICOM')
		if patientdirectory=='select':
			# Choose patient directory if not already stated
			root = Tk()
			#root.withdraw()
			direct = filedialog.askdirectory(title="select patient directory")
			root.destroy()
			print(direct)
			self.patientdirect=direct
			direct=os.path.join(direct,'DICOM')
			self.dicomdirect=direct

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
		pprint.pprint(subdirects)

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
				self.dyninfo['tres']=float(info.AcquisitionDuration/numtimepoints)
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
			self.dyninfo['tres']=float(info[0x0051,0x100a].value.split('TA ')[1])
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
			timept=temporalpos-1 # Find temporal position (begin at zero!)
			
			if temp.Manufacturer=='Philips':
				slice=int(np.floor(i/numtimepoints)) # Work out slice number (begin at zero!)\
				dynims[:,:,slice,timept]=temp.pixel_array # Read into the right part of the array
			if temp.Manufacturer=='SIEMENS':
				slice=int(i-timept*numslices)
				#print(i,slice, timept)
				dynims[:,:,slice,timept]=temp.pixel_array

		# save this file as an npy array for next time
		np.save(os.path.join(self.patientdirect,'Analysis','dynamics.npy'),dynims)
		np.save(os.path.join(self.patientdirect,'Analysis','dyninfo.npy'),self.dyninfo)
		print("dynamic image array size is "+str(dynims.shape))
		self.dynims=dynims

	def read_T1data(self,seriestag,usefolders=range(0,5,1)):
		T1Folders=glob.glob(os.path.join(self.dicomdirect,seriestag))
		if not T1Folders:
			print('Folder not found')
			return
		print('Found:')
		pprint.pprint(T1Folders)
		T1Folders=[T1Folders[x] for x in usefolders]
		print('Using these folders:' )
		pprint.pprint(T1Folders)

		# in the first directory, read in one file to check sizes etc
		T1files=glob.glob(os.path.join(T1Folders[0],'*.dcm'))
		info=dicom.read_file(T1files[0])
		im=info.pixel_array		
		self.T1info=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8'),('TIs','f8',len(T1Folders)),('N','f8')])
		self.T1info['pixelsize']=float(info.PixelSpacing[0])
		self.T1info['TR']=float(info.RepetitionTime)
		self.T1info['N']=int(info.EchoTrainLength)
		self.T1info['FlipAngle']=float(info.FlipAngle)

		if info.Manufacturer=='SIEMENS':
			print('Warning - values for N and TR may be incorrect for Siemens dicom - check:')
			print('N = ' + str(self.T1info['N']))
			print('TR = '+str(self.T1info['TR']))
			self.T1info['N']=0
			self.T1info['TR']=0

		# Make an array of the right size to hold the images
		T1ims=np.zeros(np.array([im.shape[0],im.shape[1],len(T1files),len(T1Folders)]))
		
		# Read in the rest, checking TI for each folder (don't rely on series names if Philips)
		for y in range(0,len(T1Folders)): # for each subfolder
			T1files=glob.glob(os.path.join(self.dicomdirect,T1Folders[y],'*.dcm')) # find dicom files
			temp=dicom.read_file(T1files[0]) # read first one
			if info.Manufacturer=='SIEMENS':
				self.T1info['TIs'][0][y]=float(T1Folders[y].split('TI')[1].split('DM')[0])
				ss=1
				si=0
			if info.Manufacturer=='Philips Medical Systems':
				self.T1info['TIs'][0][y]=temp[0x2001,0x101b].value # set the TI value in T1info
				#Also if Philips, make sure the scaling of the images is done
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

		print("T1 measurement image array size is "+str(T1ims.shape))
		print('Inversion times are '+str(self.T1info['TIs']))
		self.T1data=T1ims

	# Finally, a method to read all the images
	def read_ims(self,T2wseriestag,dynseriestag,T1seriestag):
		self.read_T2w(T2wseriestag)
		self.read_T1data(T1seriestag)
		self.read_dynamics(dynseriestag)

	# methods for AIF
	#####################################################
	def read_AIF_fromfile(self):
		# read existing AIF from file in patient/Analysis directory
		print('looking for AIF file')
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','AIF.npy')):
			print('AIF file not found')
			return	
		AIF=np.load(os.path.join(self.patientdirect,'Analysis','AIF.npy'))
		self.AIF=AIF

	def get_pixel_AIF(self):
		# Display dynamics and choose slices from which AIF should be extracted
		# First job - find AIF peak image by finding the index of the maximum pixel value within the first 50 frames for the central slice
		midslice=np.floor(self.dynims.shape[2]/2)
		peakframe=np.argmax(np.amax(np.amax(self.dynims[:,:,midslice,0:50],0),0))
		h2=plt.figure()
		plt.subplot(3,3,1)
		plt.imshow(self.dynims[:,:,3,peakframe])
		plt.title('3')
		plt.subplot(3,3,2)
		plt.imshow(self.dynims[:,:,4,peakframe])
		plt.title('4')
		plt.subplot(3,3,3)
		plt.imshow(self.dynims[:,:,5,peakframe])
		plt.title('5')
		plt.subplot(3,3,4)
		plt.imshow(self.dynims[:,:,6,peakframe])
		plt.title('6')
		plt.subplot(3,3,5)
		plt.imshow(self.dynims[:,:,7,peakframe])
		plt.title('7')
		plt.subplot(3,3,6)
		plt.imshow(self.dynims[:,:,8,peakframe])
		plt.title('8')
		plt.subplot(3,3,7)
		plt.imshow(self.dynims[:,:,9,peakframe])
		plt.title('9')
		plt.subplot(3,3,8)
		plt.imshow(self.dynims[:,:,10,peakframe])
		plt.title('10')
		plt.subplot(3,3,9)
		plt.imshow(self.dynims[:,:,11,peakframe])
		plt.title('11')



		ids=input("Input the id number for the slices to be used, separated by commas: ")
		ids=ids.split(',')
		ChosenSlices=[int(item) for item in ids]
		plt.close(h2)

		# Get resulting AIFs for these slices
		AIFtemp=np.zeros((2,self.dynims.shape[3],len(ChosenSlices)))
		for i in range(0,len(ChosenSlices)):
			AIFtemp[:,:,i]=getAIF.getAIF(self.dynims,ChosenSlices[i], self.dyninfo['TR']/1000, self.dyninfo['FlipAngle'])

		# Choose from plot which one to use
		h3=plt.figure()
		for ii in range(0,len(ChosenSlices)): 
			plt.plot(np.squeeze(AIFtemp[0,:,ii]),label=str(ii)+' - '+str(ChosenSlices[ii])+'L')
			plt.plot(np.squeeze(AIFtemp[1,:,ii]),label=str(ii)+' - '+str(ChosenSlices[ii])+'R')
			plt.legend()
		
		ChooseAIF=input('Which AIF should be used? Input index number, comma, 0 for L, 1 for R')
		ChooseAIF=ChooseAIF.split(',')
		AIF=AIFtemp[int(ChooseAIF[1]),:,ChooseAIF[0]]
		np.save(os.path.join(self.patientdirect,'Analysis','AIFall.npy'),AIFtemp)
		np.save(os.path.join(self.patientdirect,'Analysis','AIF.npy'),AIF) #save chosen AIF in Analysis folder as AIF.
		self.AIF=AIF
		
	# Initial processing
	#####################################################
	def make_mask(self):
		#Method to make a binary mask covering the area of interest - use T2w images
		if not hasattr(self,'T2wims'):
			print("Read in the T2w data first - patient.read_T2w()")
			return

		# Display the T2w images one at a time and mark roi if required
		numslices=self.T2wims.shape[2]
		mask=np.zeros((self.T2wims.shape),dtype='uint16')

		for jj in range(numslices):
			h3=plt.figure()
			plt.imshow(self.T2wims[:,:,jj],cmap='gray')
			useslice=input('Use this slice? 1=yes, 0=no')
			if int(useslice)==0:
				plt.close(h3)
			elif int(useslice)==1:
				roi=roipoly.roipoly()
				input('press enter to continue')
				mask[:,:,jj]=roi.getMask(self.T2wims[:,:,jj])
			else:
				print('Aborting')
				return

		self.mask=mask
		np.save(os.path.join(self.patientdirect,'Analysis','mask.npy'),mask)

	def convert_mask(self):
		#Method to convert the T2w mask to the size of the dynamic images
		if not hasattr(self,'dyninfo'):
			print('Need dynamic data to convert the mask')
			return
		
		if not hasattr(self,'mask'):
			print('No mask defined yet')
			return
		
		smallmask=np.zeros(self.dynims.shape[0:3],dtype='uint16')
		numslices=self.dynims.shape[2]

		for kk in range(numslices):
			smallmask[:,:,kk]=scipy.misc.imresize(self.mask[:,:,kk],self.dynims.shape[0:2],interp='nearest')/255

		self.dynmask=smallmask
		np.save(os.path.join(self.patientdirect,'Analysis','dynmask.npy'),smallmask)

	def load_mask(self):
		#method to load previously made dynamic mask
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','dynmask.npy')):
			print('No mask saved')
			return
		else:
			self.dynmask=np.load(os.path.join(self.patientdirect,'Analysis','dynmask.npy'))


	def make_T1map(self,seqtype='IR'): # method to do T1 fitting, seqtype is either IR or SR
		if not hasattr(self, 'T1data'):
			print('Read the T1 data first - patient.read_T1data()')
			return
		TIs=self.T1info['TIs'][0]
		TR=self.T1info['TR'][0]
		N=self.T1info['N'][0]
		flip=self.T1info['FlipAngle'][0]
		sequenceparams=(flip,np.ceil(N/2),TR,N,4000)
		
		map=np.zeros(self.dynmask.shape)
		#self.T1map=np.zeros(self.dynmask.shape)

		data=np.zeros((self.dynmask.shape+(1,)))
		data[:,:,:,0]=self.dynmask
		data=np.tile(data,(1,1,1,len(TIs)))
		data=self.T1data*data

		for sl in range(self.T1data.shape[2]): #for each slice
			print(sl)
			if np.sum(self.dynmask[:,:,sl])!=0: #if there are pixels in the slice
				for i in range(self.T1data.shape[0]):
					for j in range(self.T1data.shape[1]):
						curve=np.squeeze(data[i,j,sl,:])
						if np.sum(curve)!=0:
							curve=curve/np.amax(curve)
							if seqtype=='IR':
								fit=IRTurboFLASH.fittingfun(TIs,sequenceparams,curve)
							else:
								fit=SRTurboFLASH.fittingfun(TIs,TR,flip,np.ceil(N/2),curve)
							map[i,j,sl]=fit.x[0]
		# Save the map
		np.save(os.path.join(self.patientdirect,'Analysis','T1map.npy'),map)
		self.T1map=map

	def load_T1map(self):
		#method to load previously made dynamic mask
		if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','T1map.npy')):
			print('No T1 map saved')
			return
		else:
			self.T1map=np.load(os.path.join(self.patientdirect,'Analysis','T1map.npy'))

	
	def get_iAUC(self, baselinepts=10):
		#method to calculate iAUC for the pixels within the mask
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
								if np.isnan(np.sum(uptakeConc))==0 and np.sum(uptakeConc)>0:
									iAUC[i,j,sl]=np.trapz(uptakeConc,dx=tres)
			plt.figure()
			plt.imshow(iAUC[:,:,sl])

		# Save the map
		np.save(os.path.join(self.patientdirect,'Analysis','iAUC.npy'),iAUC)
		self.iAUC=iAUC

		def load_iAUC(self):
			#method to load previously made iAUC map
			if not os.path.isfile(os.path.join(self.patientdirect,'Analysis','iAUC.npy')):
				print('No iAUC saved')
				return
			else:
				self.iAUC=np.load(os.path.join(self.patientdirect,'Analysis','iAUC.npy'))


	# Fitting
	#####################################################
	def fit_Tofts(self):
		pass

	def fit_ExtTofts(self):
		pass

	def fit_2CXM(self,SIflag,save):
		# To fit SI, set SIflag to 1
		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF available')
			return

		if not hasattr(self,'T1map'):
			print('No T1map available')
			return

		if not hasattr(self,'dynmask'):
			print('No mask available')
			return

		self.TwoCXMfitConc=np.zeros((self.dynmask.shape+(7,)))
		self.TwoCXMfitSI=np.zeros((self.dynmask.shape+(7,)))
		
		TR=self.dyninfo['TR']/1000
		flip=self.dyninfo['FlipAngle']
		self.t=np.arange(0,self.dyninfo['tres'][0]*self.dyninfo['numtimepoints'][0],self.dyninfo['tres'][0])

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
							fit=TwoCXM.TwoCXMfittingSI(self.t, self.AIF/0.6, uptake, None, 15, TR, flip, T1base/1000,[0.03,0.3])
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
									TwoCXMfitConc=TwoCXM.TwoCXMfittingConc(self.t, self.AIF/0.6, uptakeConc, None)
									self.TwoCXMfitConc[i,j,sl,:]=TwoCXMfitConc
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','TwoCXMfitConcmaps.npy'),self.TwoCXMfitConc)

	def fit_2CUM(self, SIflag, save):
		# To fit SI, set SIflag to 1
		# To save maps, set save flag to 1
		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF available')
			return

		if not hasattr(self,'T1map'):
			print('No T1map available')
			return

		if not hasattr(self,'dynmask'):
			print('No mask available')
			return

		self.TwoCUMfitConc=np.zeros((self.dynmask.shape+(6,)))
		self.TwoCUMfitSI=np.zeros((self.dynmask.shape+(6,)))
		
		TR=self.dyninfo['TR']/1000
		flip=self.dyninfo['FlipAngle']
		self.t=np.arange(0,self.dyninfo['tres'][0]*self.dyninfo['numtimepoints'][0],self.dyninfo['tres'][0])

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
							fit=TwoCUM.TwoCUMfittingSI(self.t, self.AIF/0.6, uptake, None, 15, TR, flip, T1base/1000,[0.03,0.3])
							self.TwoCUMfitSI[i,j,sl,:]=fit
			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','TwoCUMfitSImaps.npy'),self.TwoCUMfitSI)

		else:
			for sl in range(self.dynims.shape[2]):
					for i in range(self.dynims.shape[0]):
						for j in range(0,self.dynims.shape[1]):
							if self.dynmask[i,j,sl]==1:
								uptakeConc=FLASH.SI2Conc(self.dynims[i,j,sl,:],TR,flip,self.T1map[i,j,sl]/1000,15,None)
								print(i,j,sl)
								if np.isnan(np.sum(uptakeConc))==0:
									TwoCUMfitConc=TwoCUM.TwoCUMfittingConc(self.t, self.AIF/0.6, uptakeConc, None)
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
		self.t=np.arange(0,self.dyninfo['tres'][0]*self.dyninfo['numtimepoints'][0],self.dyninfo['tres'][0])

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
								fit=AATH.AATHfittingSI(self.t, self.AIF/0.6, uptake, None, 15, TR, flip, T1base/1000)
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
									AATHfitConc=AATH.AATHfittingConc(self.t, self.AIF/0.6, uptakeConc, None)
									self.AATHfitConc[i,j,sl,:]=AATHfitConc

			if save==1:
				np.save(os.path.join(self.patientdirect,'Analysis','AATHfitConcmaps.npy'),self.AATHfitConc)


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



	def write_maps(self,dynfoldertag,mapname):
		# Method to save specified maps to dicom
		DynFolder=glob.glob(os.path.join(self.dicomdirect,dynfoldertag))
		print(os.path.join(self.dicomdirect,dynfoldertag))
		print(DynFolder)
		os.chdir(DynFolder[0])
		
		filenames=glob.glob('*.dcm')
		numfiles=len(filenames)
		maps=getattr(self,mapname)

		for i in range(maps.shape[2]):
			tmp=dicom.read_file(filenames[int(numfiles/(maps.shape[2]))*i])
			Fpimage=abs(maps[:,:,i,1]*60*100)
			Fpimage=Fpimage*(Fpimage<100)
			#print(Fpimage.dtype)
			Fpimage=Fpimage.astype('uint16')
			tmp.PixelData=Fpimage.tostring()
			dicom.write_file('Fp_'+str(i+1)+'.dcm',tmp)
			PSimage=abs(-1*np.log(1-maps[:,:,i,0])*maps[:,:,i,1])*60*100
			PSimage=PSimage*(PSimage<100)
			PSimage=PSimage.astype('uint16')
			tmp.PixelData=PSimage.tostring()
			dicom.write_file('PS_'+str(i+1)+'.dcm',tmp)
			Ktransimage=abs(maps[:,:,i,0]*maps[:,:,i,1])*60*100
			Ktransimage=Ktransimage*(Ktransimage<100)
			Ktransimage=Ktransimage.astype('uint16')
			tmp.PixelData=Ktransimage.tostring()
			dicom.write_file('Ktrans_'+str(i+1)+'.dcm',tmp)

		newfiles=glob.glob('Fp*')
		os.mkdir('Fp')
		for x in newfiles:
			os.rename(x,os.path.join('Fp',x))
		shutil.move('Fp','..')
		
		newfiles=glob.glob('PS*')
		os.mkdir('PS')
		for x in newfiles:
			os.rename(x,os.path.join('PS',x))
		shutil.move('PS','..')

		newfiles=glob.glob('Ktrans*')
		os.mkdir('Ktrans')
		for x in newfiles:
			os.rename(x,os.path.join('Ktrans',x))
		shutil.move('Ktrans','..')

	def fit_TH(self):
		pass







