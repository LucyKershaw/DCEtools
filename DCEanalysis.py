# Main file for defining the patient class to do DCE processing

# Imports
import numpy as np
import dicom
import matplotlib.pyplot as plt
import os
import glob
import PIL
import scipy.misc
import SRTurboFLASH
import scipy.optimize
import TwoCXM
import AATH
import FLASH

class patient(object): # patient inherits from the object class
	
	def __init__(self,studynumber,patientnumber,visitnumber,visitdate,hct=0.42): # set the patient number, visit number and hct for these data
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
		print('ROIs:')
		if hasattr(self,'rois'):
			print("    "+str(self.rois['roiname']))



	# Define methods to read in all the required images
	#####################################################
	def read_T2w(self): #method to read the T2w image files
		# change to dicom directory
		os.chdir(os.path.join(self.studyrootdirect,self.dicomdirect))
		# find .dcm files
		filenames=glob.glob('T2-w hr\*.dcm')
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


	def read_dynamics(self):
		# change to DynamicSeries directory
		os.chdir(os.path.join(self.studyrootdirect,self.dicomdirect,"DynamicSeries"))
		# get the timepoint directory names, sorted in the right order
		dyndirects=sorted(glob.glob('*'),key=int)
		print("Found "+str(len(dyndirects))+" time point directories")
		# in first directory, read in one slice to work out size and fill out info variables
		dynfiles=glob.glob(dyndirects[0]+"\*")
		info=dicom.read_file(dynfiles[0])
		im=info.pixel_array
		self.dyninfo=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8'),('tres','f8')])
		self.dyninfo['pixelsize']=float(info.PixelSpacing[0])
		self.dyninfo['TR']=float(info.RepetitionTime)
		self.dyninfo['FlipAngle']=float(info.FlipAngle)

		dynims=np.zeros(np.array([im.shape[0],im.shape[1],len(dynfiles),len(dyndirects)]))
		#For each directory, read files into the array
		for i in range(0,len(dyndirects)):
			dynfiles=glob.glob(dyndirects[i]+"\*")
			for j in range(0,len(dynfiles)):
				temp=dicom.read_file(dynfiles[j])
				imnum=temp.InstanceNumber
				dynims[:,:,imnum-1,i]=temp.pixel_array
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
		info=dicom.read_file(TIfiles[0])
		im=info.pixel_array		
		self.T1info=np.zeros(1,dtype=[('pixelsize','f8'),('TR','f8'),('FlipAngle','f8'),('TIs','f8',len(TIdirects)),('n','f8')])
		self.T1info['pixelsize']=float(info.PixelSpacing[0])
		#self.T1info['TR']=float(info.RepetitionTime)
		self.T1info['TR']=float(4) #Echo spacing is 4 ms in HN data - check!
		self.T1info['n']=float(5) #HN data - check! 
		self.T1info['FlipAngle']=float(info.FlipAngle)
		self.T1info['TIs']=np.array([float(item.split('=')[1]) for item in TIdirects])

		SRims=np.zeros(np.array([im.shape[0],im.shape[1],len(TIfiles),len(TIdirects)]))
		# for each TI directory, read files into the array
		for i in range(0,len(TIdirects)):
			TIfiles=glob.glob(TIdirects[i]+"\*")
			for j in range(0,len(TIfiles)):
				temp=dicom.read_file(TIfiles[j])
				imnum=temp.InstanceNumber
				SRims[:,:,imnum-1,i]=temp.pixel_array
		print("T1 measurement image array size is "+str(SRims.shape))
		self.T1data=SRims

	# Finally, a method to read all the images
	def read_ims(self):
		self.read_T2w()
		self.read_T1data()
		self.read_dynamics()

	# Define methods to read in and convert the ROI files, and also to show which ROIs have been read
	#####################################################
	def read_rois(self,readall=1):
		# readall is flag to read all .roi files in the folder.  if you want to choose, set readall to 0 instead
		# Check for T2w images - need their size to get the roi array size
		if not hasattr(self,'T2wims'):
			print("Please read the T2w images first")
			return
		if not hasattr(self,'dynims'):
			print("Please read the dynamic images first")
			return

		# change to roi directory and display the available ROI files in the ROI directory
		os.chdir(os.path.join(self.studyrootdirect,self.roidirect))
		roifiles=glob.glob('*.roi')
		print("ROI files found:")
		print(roifiles)

		# deal with choosing files if requested
		if readall==0:
			[print(item) for item in enumerate(roifiles)]
			ids=input("Input the id number for the chosen file(s), separated by commas: ")
			ids=ids.split(',')
			roifiles=[roifiles[int(item)] for item in ids]
			print("Chosen rois:")
			[print(item) for item in roifiles]

		#make the appropriate structured array to put the read rois into
		self.rois=np.zeros(len(roifiles),dtype=[('roiname','a60'),('hiresarray','u2',self.T2wims.shape),('dynresarray','u2',self.dynims.shape[0:3]),('T1','f8')])

		# go through each file, reading and converting to arrays

		for i in range (0,len(roifiles)):
			#Read the file  - unsigned short little endian
			roi=np.fromfile(roifiles[i],np.uint16)
			# Set flag for extended format
			extendedformat=0
			# Make a results array to hold the roi
			mask=np.zeros(np.array(self.T2wims.shape))
			smallmask=np.zeros(np.array(self.dynims.shape[0:3]))

			# Now begin parsing the array - the structure of the .roi file from MRIcro is [Slice index, words in slice, 
			# begin index, run length, begin index, run length..... slice index...]
			start_pos=0
			end_pos=1

			while end_pos < roi.size:
				SliceIndex=roi[start_pos]
				# Check for extended format using slice index
				if (SliceIndex-(256*256)/2)>0:
					extendedformat=1
					SliceIndex=SliceIndex-(256*256)/2
				# Assign temporary array to put the values into
				slice_tmp=np.zeros(np.array(self.T2wims.shape[0:2]))
				WordsInSlice=roi[start_pos+1] # Get the number of words - that's in the next position
				end_pos=start_pos+WordsInSlice # Set the place to stop so that only words for this slice are read
				tmp=roi[start_pos+2:end_pos] # read all the words for this slice (+2 to avoid reading sliceindex and wordsinslice again)
				Indexes=np.zeros([2,tmp.size/2],dtype=np.uint16) # Clear these from last time round the loop
				longIndexes=np.zeros([2,tmp.size/2],dtype=np.uint64)
				
				Indexes[0,:]=tmp[0::2] # first row contains the begin indices (start, stop, step)
				Indexes[1,:]=tmp[1::2] # second row contains the run lengths
				start_pos=start_pos+WordsInSlice # Set the new starting position at the next SliceIndex
				if extendedformat: # This is slightly unclear, but extended format involves a bit of shifting of bits to make it work
					addit=np.right_shift(Indexes[1,:],12)*(256*256) # shift right 12 bit, mult. 65536
					subt=np.right_shift(Indexes[1,:],12)*(256*16) # same, but mult. by 4096
					oldInd=Indexes;
					longIndexes[0,:]=Indexes[0,:]+addit 	# add to begin index
					longIndexes[1,:]=Indexes[1,:]-subt # subtract offset * 4096
					Indexes=np.uint64(Indexes)
					Indexes=longIndexes
				for k in range(0,np.shape(Indexes)[1]):
					slice_tmp.ravel()[Indexes[0,k]-1:(Indexes[0,k]-1+Indexes[1,k])]=1 # set mask voxels to 1 using linear indexing, -1 to take account of starting at 0
				#plt.imshow(np.flipud(slice_tmp),cmap='gray',vmin=0,vmax=1)
				# finally, put this into the full mask
				mask[:,:,SliceIndex-1]=np.flipud(slice_tmp) # flip up-down to make this the same orientation as a read-in dicom image
				smallmask[:,:,SliceIndex-1]=scipy.misc.imresize(np.flipud(slice_tmp),self.dynims.shape[0:2],interp='nearest')/255

			# put the mask into the rois structured array
			self.rois['roiname'][i]=roifiles[i]
			self.rois['hiresarray'][i]=mask
			self.rois['dynresarray'][i]=smallmask

	# methods for AIF
	#####################################################
	def read_AIF_fromfittingfile(self,AIFdirectory='E:\HeadAndNeck\HeadNeckNewPts\Fitting\Tumour'):
		# read existing AIF from fitting file in given directory
		# Usually called P50_pre.txt, etc, found in third column in conc units assuming hct of 0.42
		# and r1 of 4.5

		if self.visitnumber==1:
			visittext='pre'
		if self.visitnumber==2:
			visittext='post'

		filename='P'+str(self.studynumber)+'_'+visittext+'.txt'
		print('looked for a file called '+filename)
		AIFfile=np.loadtxt(os.path.join(AIFdirectory,filename))
		# Convert back to delta R1 and use patient hct
		self.AIF=(AIFfile[:,2]*4.5)*(1-0.42)/(1-self.hct)
		self.t=AIFfile[:,1] # time in minutes

	def getpixelAIF(self):
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


	# Display methods
	#####################################################
	#def 
