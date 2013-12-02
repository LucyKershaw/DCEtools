# Main file for defining the patient class to do DCE processing

# Imports
import numpy as np
import dicom
import matplotlib.pyplot as plt
import os
import glob

class patient(object): # patient inherits from the object class
	
	def __init__(self,studynumber,patientnumber,visitnumber,visitdate,hct): # set the patient number, visit number and hct for these data
		self.studyrootdirect=os.path.normpath("E:/HeadAndNeck/HeadNeckNewPts/Patients/") # this will be head and neck data
		self.hct=hct
		self.studynumber=studynumber
		self.patientnumber=patientnumber
		self.visitnumber=visitnumber
		self.visitdate=visitdate
		self.dicomdirect=os.path.normpath("Patient "+str(studynumber)+"/"+str(visitdate)+"/Dicom")# Construct the expected name for the dicom directory
		self.roidirect=os.path.normpath("Patient "+str(studynumber)+"/"+str(visitdate)+"/ROIs")# Construct the expected name for the roi directory
		print("Patient "+str(studynumber)+" visit "+str(visitnumber)+" created.  Imaging data expected in:")
		print(os.path.join(self.studyrootdirect,self.dicomdirect)+" Does this exist? "+str(os.path.isdir(os.path.join(self.studyrootdirect,self.dicomdirect)))) # Check these exist
		print(os.path.join(self.studyrootdirect,self.roidirect)+" Does this exist? "+str(os.path.isdir(os.path.join(self.studyrootdirect,self.roidirect))))

	# Define methods to read in all the required images
	#####################################################
	def read_T2w(self): #method to read the T2w image files
		# change to dicom directory
		os.chdir(os.path.join(self.studyrootdirect,self.dicomdirect))
		# find .dcm files
		filenames=glob.glob('T2-w hr\*.dcm')
		numfiles=len(filenames)
		print("Reading "+str(numfiles)+" files")
		# read the first file to find out size
		im=dicom.read_file(filenames[0]).pixel_array

		T2wims=np.zeros(np.array([im.shape[0],im.shape[1],numfiles]))
		for i in range(0,numfiles-1):
			T2wims[:,:,i]=dicom.read_file(filenames[i]).pixel_array
		self.T2wims=T2wims


	def read_dynamics(self):
		# change to DynamicSeries directory
		os.chdir(os.path.join(self.studyrootdirect,self.dicomdirect,"DynamicSeries"))
		# get the timepoint directory names, sorted in the right order
		dyndirects=sorted(glob.glob('*'),key=int)
		print("Found "+str(len(dyndirects))+" time point directories")
		# in first directory, read in one slice to work out size
		dynfiles=glob.glob(dyndirects[0]+"\*")
		im=dicom.read_file(dynfiles[0]).pixel_array
		dynims=np.zeros(np.array([im.shape[0],im.shape[1],len(dynfiles),len(dyndirects)]))
		#For each directory, read files into the array
		for i in range(0,len(dyndirects)):
			dynfiles=glob.glob(dyndirects[i]+"\*")
			for j in range(0,len(dynfiles)):
				dynims[:,:,j,i]=dicom.read_file(dynfiles[j]).pixel_array
		print("dynamic image array size is "+str(dynims.shape))
		self.dynims=dynims

	def read_T1data(self):
		# change to SR directory
		os.chdir(os.path.join(self.studyrootdirect,self.dicomdirect,"SR"))
		# get the TI directory names, sorted correctly
		TIdirects=glob.glob('*')
		TIdirects.sort(key=lambda x: int(x[3:]))
		print("Found these TI directories - ")
		print(TIdirects)
		# in the first directory, read in one file to check sizes
		TIfiles=glob.glob(TIdirects[0]+"\*")
		im=dicom.read_file(TIfiles[0]).pixel_array
		SRims=np.zeros(np.array([im.shape[0],im.shape[1],len(TIfiles),len(TIdirects)]))
		# for each TI directory, read files into the array
		for i in range(0,len(TIdirects)):
			TIfiles=glob.glob(TIdirects[i]+"\*")
			for j in range(0,len(TIfiles)):
				SRims[:,:,j,i]=dicom.read_file(TIfiles[j]).pixel_array
		print("T1 measurement image array size is "+str(SRims.shape))
		self.T1data=SRims

	# Define methods to read in and convert the ROI files
	#####################################################
	def read_ROI(self):
		# display the available ROI files
		

	# methods for AIF
	#####################################################
	def read_AIF(self):
		pass # read existing AIF from file

	def getpixelAIF(self):
		pass # get an AIF from the dynamic data
		
	# Initial processing
	#####################################################
	def applyROI(self):
		pass # method to extract data about the ROI alone - resample the ROI to dynamic matrix and get T1 and dynamic SI curves

	def fit_T1(self): # method to do the T1 fitting
		pass # method to fit T1 for ROI

	def SIconvert(self): 
		pass # convert the dynamic curve to concentration

	# Fitting
	#####################################################
	def fit_2CXM(self):
		pass # Fit the 2CXM to the dynamic curve

	# Display methods
	#####################################################
	#def 
