
# Module to resample one image on the grid of another

import glob
import numpy as np
import dicom
import os
import scipy.interpolate
import scipy.ndimage.morphology


# Resample one dicom (orig_image) into the space of another (target_image)
# Inputs are the image folders containing all the slices (need to know position of first and last slice)
# Optionally include an already-read nparray for the original image
# CURRENTLY ASSUMES AXIAL TO AXIAL

def resample_dicom(target_image_folder,orig_image_folder,orig_nparray=np.zeros((3,3)),target_slices=0,orig_slices=0): # Resample one dicom (orig_image) into the space of another (target_image), optionally include an already-read nparray for the original image

	# Read necessary info from dicom files
	# Read the IPP, pixel spacing and matrix from the dicom header for the two image sets
	target_files=glob.glob(os.path.join(target_image_folder,'*.dcm'))
	orig_files=glob.glob(os.path.join(orig_image_folder,'*.dcm'))

	#ImagePositionPatient
	target_IPP0=np.array(dicom.read_file(target_files[0]).ImagePositionPatient)
	target_IPPN=np.array(dicom.read_file(target_files[-1]).ImagePositionPatient)
	orig_IPP0=np.array(dicom.read_file(orig_files[0]).ImagePositionPatient)
	orig_IPPN=np.array(dicom.read_file(orig_files[-1]).ImagePositionPatient)

	#ImageOrientationPatient
	target_IOP=np.array(dicom.read_file(target_files[0]).ImageOrientationPatient)
	orig_IOP=np.array(dicom.read_file(orig_files[0]).ImageOrientationPatient)
	
	#PixelSpacing
	target_PixelSpacing=dicom.read_file(target_files[0]).PixelSpacing
	orig_PixelSpacing=dicom.read_file(orig_files[0]).PixelSpacing

	#Num rows, cols, slices
	target_rows=dicom.read_file(target_files[0]).Rows
	target_cols=dicom.read_file(target_files[0]).Columns
	if target_slices==0:
		target_slices=int(len(target_files)/dicom.read_file(target_files[-1]).AcquisitionNumber)
	print('Target image array')
	print(str(target_rows)+','+str(target_cols)+','+str(target_slices))
	orig_rows=dicom.read_file(orig_files[0]).Rows
	orig_cols=dicom.read_file(orig_files[0]).Columns
	if orig_slices==0:
		orig_slices=int(len(orig_files)/dicom.read_file(orig_files[-1]).AcquisitionNumber)
	print('Original image array')
	print(str(orig_rows)+','+str(orig_cols)+','+str(orig_slices))


	#Make affines
	target_A=np.zeros((4,4))
	target_A[0:3,0]=target_IOP[3:6]*target_PixelSpacing[0]
	target_A[0:3,1]=target_IOP[0:3]*target_PixelSpacing[1]
	target_A[0:3,3]=target_IPP0
	target_A[3,3]=1
	target_A[0:3,2]=(target_IPP0-target_IPPN)/(1-target_slices)


	orig_A=np.zeros((4,4))
	orig_A[0:3,0]=orig_IOP[3:6]*orig_PixelSpacing[0]
	orig_A[0:3,1]=orig_IOP[0:3]*orig_PixelSpacing[1]
	orig_A[0:3,3]=orig_IPP0
	orig_A[3,3]=1
	orig_A[0:3,2]=(orig_IPP0-orig_IPPN)/(1-orig_slices)		


	# Now need to make grids from the affines
	target_r=np.zeros((target_rows,target_cols,target_slices))
	target_c=np.zeros((target_rows,target_cols,target_slices))
	target_s=np.zeros((target_rows,target_cols,target_slices))
	#print(MPRAGEr.shape)

	for r in range(target_rows):
		for c in range(target_cols):
			for s in range(target_slices):
				coords=target_A.dot([r,c,s,1])
				target_r[r,c,s]=coords[0]
				target_c[r,c,s]=coords[1]
				target_s[r,c,s]=coords[2]

	print(target_r.dtype)

	orig_r=np.zeros((orig_rows,orig_cols,orig_slices))
	orig_c=np.zeros((orig_rows,orig_cols,orig_slices))
	orig_s=np.zeros((orig_rows,orig_cols,orig_slices))
	#print(Dynr.shape)

	for r in range(orig_rows):
		for c in range(orig_cols):
			for s in range(orig_slices):
				coords=orig_A.dot([r,c,s,1])
				orig_r[r,c,s]=coords[0]
				orig_c[r,c,s]=coords[1]
				orig_s[r,c,s]=coords[2]

	print(orig_r.dtype)	

	if np.sum(orig_nparray)==0: #If nparray isn't provided, read in the original images (ASSUMING NOT DYNAMIC)
		print("Reading original images from file")
		orig_nparray=np.zeros(np.array([orig_rows,orig_cols,orig_slices]))
		for i in range(0,orig_slices):
			temp=dicom.read_file(orig_files[i])
			imnum=temp.InstanceNumber
			orig_nparray[:,:,imnum-1]=temp.pixel_array

	# Resample the data

	print('Resampling...')
	# Interpolate original image data onto target grid using nearest neighbour interpolation
	resampled_image=scipy.interpolate.griddata((orig_r.flatten(),orig_c.flatten(),orig_s.flatten()),orig_nparray.flatten(),(target_r,target_c,target_s),method='nearest')

	return(resampled_image)

	