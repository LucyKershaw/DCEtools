# Main file for defining the patient class to do DCE processing

class patient(object): # patient inherits from the object class
	# Imports
	import numpy as np
	import dicom

	def __init__(self,patientnumber,visitnumber,hct): # set the patient number, visit number and hct for these data
		self.patientnumber=patientnumber
		self.visitnumber=visitnumber
	
	# Define methods to read in all the required images
	#####################################################
	def read_T2w(self):
		pass # add code to read T2w images
	def read_dynamics(self):
		pass # add code to read dynamics
	def read_T1data(self):
		pass # add code to read T1 data

	# Define methods to read in and convert the ROI files
	#####################################################
	def read_ROI(self):
		pass

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
	def 
