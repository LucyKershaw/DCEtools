# Function to convert mricro .roi files to binary masks

def roi2mask(filename,x,y,z):
	# filedims is the size of the images that the ROI was marked on
	import numpy as np
	import pylab as plt
	#Read the file  - unsigned short little endian
	roi=np.fromfile(filename,np.uint16)
	# Set flag for extended format
	extendedformat=0
	# Make a results array to hold the roi
	mask=np.zeros([x,y,z])

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
		slice_tmp=np.zeros([x,y])
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
		print(Indexes)
		for i in range(0,np.shape(Indexes)[1]):
			slice_tmp.ravel()[Indexes[0,i]-1:(Indexes[0,i]-1+Indexes[1,i])]=1 # set mask voxels to 1 using linear indexing, -1 to take account of starting at 0
		plt.imshow(np.flipud(slice_tmp),cmap='gray',vmin=0,vmax=1)
		# finally, put this into the full mask
		mask[:,:,SliceIndex-1]=np.flipud(slice_tmp) # flip up-down to make this the same orientation as a read-in dicom image
	return mask # return the mask as a numpy array 

