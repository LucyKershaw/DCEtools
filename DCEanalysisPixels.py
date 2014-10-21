# Main file for defining the patient class to do DCE processing pixel by pixel

# Imports
import numpy as np
import dicom
import matplotlib.pyplot as plt
import os
import glob
from tkinter import Tk
from tkinter import filedialog

import scipy.misc
import SRTurboFLASH
import IRTurboFLASH
import scipy.optimize
import TwoCXM
import AATH
import FLASH
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
		for i in range(0,numfiles-1):
			temp=dicom.read_file(filenames[i])
			imnum=temp.InstanceNumber
			T2wims[:,:,imnum-1]=temp.pixel_array
		self.T2wims=T2wims


	def read_dynamics(self,seriestag):
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
		# Check for Philips data - file numbering different for Siemens
		if info.Manufacturer!='Philips Medical Systems':
			print('This is not Philips dynamic data - need to rewrite for Siemens!')
			return

		# Find pixel sizes and other required dynamic info
		self.dyninfo=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8'),('tres','f8'),('numtimepoints','i4'),('numslices','i4')])
		self.dyninfo['pixelsize']=float(info.PixelSpacing[0])
		self.dyninfo['TR']=float(info.RepetitionTime)
		self.dyninfo['FlipAngle']=float(info.FlipAngle)
		numtimepoints=int(info.TemporalPositionIdentifier)
		self.dyninfo['numtimepoints']=numtimepoints
		self.dyninfo['tres']=float(info.AcquisitionDuration/numtimepoints)
		numslices=int(len(dynfiles)/numtimepoints)
		self.dyninfo['numslices']=numslices

		# Make an array to hold the dynamic data
		dynims=np.zeros(np.array([im.shape[0],im.shape[1],numslices,numtimepoints]))
		
		# Read files into the array
		for i in range(0,len(dynfiles)):
			temp=dicom.read_file(dynfiles[i]) # Read file
			timept=temp.TemporalPositionIdentifier-1 # Find temporal position (begin at zero!)
			slice=int(np.floor(i/numtimepoints)) # Work out slice number (begin at zero!)
			dynims[:,:,slice,timept]=temp.pixel_array # Read into the right part of the array

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
		self.T1info=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8'),('TIs','f8',len(T1Folders)),('n','f8')])
		self.T1info['pixelsize']=float(info.PixelSpacing[0])
		self.T1info['TR']=float(info.RepetitionTime)
		self.T1info['n']=int(info.EchoTrainLength/2)
		self.T1info['FlipAngle']=float(info.FlipAngle)

		if info.Manufacturer=='SIEMENS':
			print('Warning - values for n and TR are probably incorrect for Siemens dicom')
			print('Setting these to zero, set them manually')
			self.T1info['n']=0
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
		if not os.path.isfile(os.path.join(self.patientdirect,Analysis,'AIF.txt')):
			print('AIF file not found')
			return	
		AIFfile=np.genfromtxt(os.path.join(self.patientdirect,Analysis,'AIF.txt'))
		self.AIF=AIFfile

	def get_pixel_AIF(self):
		pass # get an AIF from the dynamic data
		
	# Initial processing
	#####################################################
	def get_T1curves(self):
		# method to extract T1 curve from T1 images for loaded rois
		if not hasattr(self,'rois'):
			print("Read the roi files first - patient.read_rois()")
			return
		if not hasattr(self,'T1data'):
			print("Read in the T1 data first - patient.read_T1data")
			return
		curves=np.zeros((self.T1data.shape[3],len(self.rois)))
		for i in range(0,len(self.rois)):
			for j in range(0,self.T1data.shape[3]):
				curves[j,i]=np.sum(self.T1data[:,:,:,j]*self.rois['dynresarray'][i])/np.sum(self.rois['dynresarray'][i])
		self.T1curves=curves

	def get_SIcurves(self):
		# method to extract SI curves for the loaded rois
		if not hasattr(self,'rois'):
			print("Read the roi files first - patient.read_rois()")
			return
		if not hasattr(self,'dynims'):
			print("Read in the dynamic data first - patient.read_dynamics")
			return
		curves=np.zeros((self.dynims.shape[3],len(self.rois)))
		

		for i in range(0,len(self.rois)): #for each roi..
			# add a fourth dimension
			mask=self.rois['dynresarray'][i]
			mask.shape=mask.shape+(1,) 
			# duplicate this mask for all time points using tile
			ntimepoints=self.dynims.shape[3]
			bigmask=np.tile(mask,(1,1,1,ntimepoints))
			#multiply by whole dynamic array, sum over timepoints and divide by mask sum
			curves[:,i]=np.sum((bigmask*self.dynims),(0,1,2))/np.sum(mask)

		self.SIcurves=curves

	def fit_T1s(self,plotfit=0): # method to do the T1 fitting, set plotfit to 1 if plotting required
		if not hasattr(self, 'T1curves'):
			print('Extract the T1 curves first - patient.get_T1curves()')
			return
		TIs=self.T1info['TIs'][0]
		TR=self.T1info['TR']
		n=self.T1info['n']
		flip=self.T1info['FlipAngle']

		for i in range(0,self.T1curves.shape[1]):
			data=self.T1curves[:,i]
			fit=SRTurboFLASH.fittingfun(TIs,TR,flip,n,data)
			if plotfit==1:
				plt.plot(TIs,data,'x')
				plt.plot(TIs,SRTurboFLASH.SIeqn(fit.x,TIs,TR,flip,n))
			self.rois['T1'][i]=fit.x[0]

	def SIconvert(self, baselinepts=10): 
		#Check we have rois and T1s
		if not hasattr(self,'rois'):
			print('Read in rois first - patient.read_rois()')
			return
		if sum(self.rois['T1'])==0:
			print('fit T1 first - patient.fit_T1s()')
			return
		# Convert flip angle to radians
		rflip=self.dyninfo['FlipAngle']*np.pi/180
		#extract TR
		TR=self.dyninfo['TR']/1000
		Conccurves=np.zeros(self.SIcurves.shape)

		for i in range(0,len(self.rois)):
			SIcurve=self.SIcurves[:,i]
			T1base=self.rois['T1'][i]
			# Convert T1 to R1 in s^-1
			R1base=1/(T1base/1000)
			print(R1base)
			# extract baseline SI and calculate M0
			base=np.mean(SIcurve[0:baselinepts])
			print(base)
			M0=base*(1-np.cos(rflip)*np.exp(-1*TR*R1base))/(np.sin(rflip)*(1-np.exp(-1*TR*R1base)))
			print(M0)
			# Now calculate the R1 curve
			R1=np.log(((M0*np.sin(rflip))-SIcurve)/(M0*np.sin(rflip)-(SIcurve*np.cos(rflip))))*(-1/TR)
			# And finally the delta R1 curve
			Conccurves[:,i]=R1-R1base

		self.Conccurves=Conccurves

	def preprocess(self):
		# function to extract and fit T1 curves, extract SI curves and convert
		self.get_T1curves()
		self.get_SIcurves()
		self.fit_T1s()
		self.SIconvert()

	# Fitting
	#####################################################
	def fit_Tofts(self):
		pass

	def fit_ExtTofts(self):
		pass

	def fit_2CXM(self,SIflag):
		# To fit SI, set SIflag to 1
		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF - if reading from file, patient.read_AIF_fromfittingfile')
			return
		TwoCXMfitConc=np.zeros([6,len(self.rois)])
		TwoCXMfitSI=np.zeros([6,len(self.rois)])
		
		if SIflag==1:
			for i in range(0,len(self.rois)):
				uptake=self.SIcurves[:,i]
				TR=self.dyninfo['TR']/1000
				flip=self.dyninfo['FlipAngle']
				T1base=self.rois['T1'][i]
				TwoCXMfitSI[:,i]=TwoCXM.TwoCXMfittingSI(self.t, self.AIF, uptake, None, 5, TR, flip, T1base/1000)
				self.TwoCXMfitSI=TwoCXMfitSI

		else:
			for i in range(0,len(self.rois)):
				uptake=self.Conccurves[:,i]
				TwoCXMfitConc[:,i]=TwoCXM.TwoCXMfittingConc(self.t,self.AIF,uptake,None) # Fit the 2CXM to the dynamic curve
				self.TwoCXMfitConc=TwoCXMfitConc


	def fit_AATH(self,SIflag):
		# To fit SI, set SIflag to 1
		# check for an AIF
		if not hasattr(self,'AIF'):
			print('No AIF - if reading from file, patient.read_AIF_fromfittingfile')
			return
		AATHfitConc=np.zeros([6,len(self.rois)])
		AATHfitSI=np.zeros([6,len(self.rois)])
		
		if SIflag==1:
			for i in range(0,len(self.rois)):
				uptake=self.SIcurves[:,i]
				TR=self.dyninfo['TR']/1000
				flip=self.dyninfo['FlipAngle']
				T1base=self.rois['T1'][i]
				AATHfitSI[:,i]=AATH.AATHfittingSI(self.t, self.AIF, uptake, None, 5, TR, flip, T1base/1000)
				self.AATHfitSI=AATHfitSI

		else:
			for i in range(0,len(self.rois)):
				uptake=self.Conccurves[:,i]
				AATHfitConc[:,i]=AATH.AATHfittingConc(self.t,self.AIF,uptake,None) # Fit the 2CXM to the dynamic curve
				self.AATHfitConc=AATHfitConc


	def fit_TH(self):
		pass


	# # Export methods
	# #####################################################
	# def export_fits(self,exportfilename):
	# 	# Remember to add the full path for the results
		



	# 	# Create the file if not already there
	# 	if not os.isfile(exportfilename):
	# 		f=open(exportfilename,'x')
	# 		f.close
	# 	else:
	# 		with open(exportfilename,'a') as csvfile:
 #    		writefile=csv.writer(csvfile)
 #    		writefile.writerow(squeeze(P42.AATHfitSI))





