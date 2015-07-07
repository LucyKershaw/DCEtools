# Function to extract contours from a structure set file

def getContours():
	import dicom
	import numpy as np
	from tkinter import Tk
	from tkinter import filedialog
	import matplotlib.pyplot as plt
	
	# Choose contour file
	root = Tk()
	#root.withdraw()
	ssfilename = filedialog.askopenfilename(title="select structure set file")
	root.destroy()
	print(ssfilename)

	#Read the structure set file
	ss=dicom.read_file(ssfilename)

	#Find number of contours in here
	numcontours=len(ss.ROIContours)

	#Print contour names and get associated numbers
	print('These contours are named in the file')	
	ROINames=[ss.StructureSetROIs[x].ROIName for x in range(0,numcontours)]
	ROINums=[ss.StructureSetROIs[y].ROINumber for y in range(0,numcontours)]
	ROIContournums=[ss.ROIContours[z].ReferencedROINumber for z in range(0,numcontours)]

	print(ROINames)
	print(ROINums) # must be matched with ss.ROIContours[i].ReferencedROINumber because there is potential for these to be stored out of order
	print(ROIContournums)

	#Now read the images that the contours correspond to
	root = Tk()
	#root.withdraw()
	imfilenames = filedialog.askopenfilenames(title="select matching image files")
	root.destroy()
	#Sort the filenames according to slice location
	IPP=[dicom.read_file(x).ImagePositionPatient for x in imfilenames]
	IPPz=[f[2] for f in IPP]
	together=zip(IPPz,imfilenames,IPP)
	together=sorted(together)
	sortedimfilenames=[y[1] for y in together]
	sortedUIDfromimfilenames=[x.split('/')[-1][3:] for x in sortedimfilenames]

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

	#Choose contour for display
	[print(item) for item in enumerate(ROINames)]
	id=input("Input the id number for the chosen contour")
	ROIchoice=ROINames[int(id)]
	print("Chosen roi:")
	print(ROIchoice)

	#Check that ROIContournum matches ROINum
	if ROINums[int(id)]!=ROIContournums[int(id)]:
		print('Warning - contours stored out of order, check this is the right structure')

	#Find the imageUIDs for this contour
	Contour=ss.ROIContours[int(id)].ContourSequence
	ContourRefdUID=[Contour[x].ContourImageSequence[0].ReferencedSOPInstanceUID for x in range(0,len(Contour))]

	#Loop through contour slices and display
	for i in range(0,len(Contour)):
		#Find index of correct slice in image stack
		ind=sortedUIDfromimfilenames.index(ContourRefdUID[i])
		#Display using correct coordinates
		plt.figure()
		plt.imshow(images[:,:,ind],extent=[sortedIPP[ind][0],sortedIPP[ind][0]+cols*PS[0],sortedIPP[ind][1]+rows*PS[0],sortedIPP[ind][1]],interpolation='nearest',cmap='gray')
		SliceContourData=Contour[i].ContourData
		plt.plot(SliceContourData[0::3],SliceContourData[1::3],'r')