# Function to rename and sort dicom files from Philips scanner


def Filesort_Philips(direct):

	import dicom
	import glob
	import os
	import shutil
	from tkinter import Tk
	from tkinter import filedialog

	# First choose a directory for sorting, if not specified
	if direct=='':
		root = Tk()
		root.withdraw()
		direct = filedialog.askdirectory(title="select patient directory for sorting")
	
	print(direct)
	direct=os.path.join(direct,'DICOM')

	# Now see what subdirectories are inside the DICOM folder in the directory
	subdirects=[x for x in os.listdir(direct) if os.path.isdir(os.path.join(direct, x))]

	# Work through the main DICOM directory and these subdirectories, reading the DICOM header
	# and changing the filenames accordingly

	# Main directory
	filenames=glob.glob(os.path.join(direct,'I*'))
	filenames.sort()
	print('Processing main directory')

	for N in filenames:
		dicom_renamePhilips(N,direct)
		
	# Subdirectories
	for M in subdirects:
		filenames=glob.glob(os.path.join(direct,M,'I*')) # find image files
		filenames.sort()
		print(len(filenames))
		print('Processing subdirectory '+M)
		for N in filenames:
			dicom_renamePhilips(N,direct) # rename them all and move to DICOM folder
		# Delete remaining files and the subdirectory itself	
		shutil.rmtree(os.path.join(direct,M),ignore_errors=True)

	# Remove files that aren't dicom
	for P in glob.glob(os.path.join(direct,'XX*')):
		os.remove(P)
	for T in glob.glob(os.path.join(direct,'PS*')):
		os.remove(T)

	# Sort into folders for each series
	sort_series(direct)

	# Add an Analysis folder in the main patient directory
	os.mkdir(os.path.join(os.path.split(direct)[0],'Analysis'))

def Filesort_Siemens(direct,collapse_dynamics=1):
	import dicom
	import glob
	import os
	import shutil
	from tkinter import Tk
	from tkinter import filedialog
	from itertools import chain
	from itertools import groupby
	from operator import itemgetter

	print(direct)
	# First choose a directory for sorting, if not specified
	if direct=='':
		root = Tk()
		root.withdraw()
		direct = filedialog.askdirectory(title="select patient directory for sorting")
	
	print(direct)
	direct=os.path.join(direct,'DICOM')

	# Main directory
	filenames=glob.glob(os.path.join(direct,'*.dcm'))
	filenames.sort()
	if len(filenames)==0:
		filenames=glob.glob(os.path.join(direct,'*.IMA'))
	print('Processing directory')

	for N in filenames:
		dicom_renameSiemens(N,direct)

	# Remove files that aren't dicom
	for P in glob.glob(os.path.join(direct,'*.txt')):
		os.remove(P)

	#If vibe images are all saved as separate series, collapse into one for all dynamic run
	if collapse_dynamics==1:
		#get a list of VIBE images
		vibefiles=glob.glob(os.path.join(direct,'*vibe*'))
		vibefiles.sort()
		vibefiles=[os.path.basename(Y) for Y in vibefiles]
		#find the seriesnumbers
		allvibenums=[int(X.split('_')[0]) for X in vibefiles]
		vibenums=set(allvibenums)    #get unique seriesnums
		vibenums=list(vibenums)   #change back to list
		#group according to consecutive series numbers
		temp=[list(map(itemgetter(1),g)) for k,g in groupby(enumerate(vibenums),lambda x:x[0]-x[1])]
		#merge lists for consecutive series groups greater than 20 (note that series 99 is commonly left out!)
		merged=list(chain.from_iterable([X for X in temp if len(X)>20]))
		if len(merged)==0: #if there are no long sets of continuous vibe seriesnumbers, sort the rest and exit
			sort_series(direct)
			return

		#Find files in the original list that have these numbers
		dynfiles=[s for s in vibefiles if int(s.split('_')[0]) in merged]
		#put into a dynamic folder with the series number from the first file
		os.mkdir(os.path.join(direct,str(min(merged))+'_'+dynfiles[0].split('_')[1])) #make folder with approriate name
		for Z in dynfiles:
			sernum=int(Z.split('_')[0])
			sernum='%03d' % sernum
			newname=sernum+'_'+Z.split('_')[1]+'_'+Z.split('_')[2]
			#print(newname)
			os.rename(os.path.join(direct,Z),os.path.join(direct,str(min(merged))+'_'+dynfiles[0].split('_')[1],newname))


	# Sort into folders for each series
	sort_series(direct)

##############################################################################################
# Function to rename Siemens dicom within series folders so that images are in sequential order

def rename_in_folders(direct, UMCUPhilips=0, rename_folders=0, UoMGE=0):
	import dicom
	import os
	import glob

	#Find existing series folders
	foldernames=glob.glob(os.path.join(direct,'*'))
	foldernames.sort()
	#print(foldernames)

	if UMCUPhilips==1:
		#For each folder...
		#find existing filenames
		for j in range(len(foldernames)):
			filenames=glob.glob(os.path.join(foldernames[j],'I*'))
			filenames.sort()
			# Then rename each file
			for i in range(len(filenames)):
				dicom_renamePhilips(filenames[i],foldernames[j])
				if i%10==0:
					print(i)
		if rename_folders==1:
			renamefolders(foldernames)
		return



	#For each folder...
	#find existing filenames
	for j in range(len(foldernames)):
		filenames=glob.glob(os.path.join(foldernames[j],'*.dcm'))
		filenames.sort()
		if UoMGE==1: # From XNAT, series folders have an extra directory inside them
			filenames=glob.glob(os.path.join(foldernames[j],'DICOM/*.dcm'))
			filenames.sort()
		# Then rename each file
		for i in range(len(filenames)):
			dicom_renameSiemens(filenames[i],foldernames[j])
			if i%10==0:
				print(i)
		if UoMGE==1:
			os.rmdir(os.path.join(foldernames[j],'DICOM'))

	if rename_folders==1:
		renamefolders(foldernames)




###############################################################################################
# Helper function to read in a header, get the info for a new filename and rename the file
# The file is then copied to the destination folder dstfolder

def dicom_renamePhilips(filename, dstfolder):
	import dicom
	import os

	im=dicom.read_file(filename) # Read the file
	if im.SOPClassUID=='1.2.840.10008.5.1.4.1.1.66':
		print('Raw data file found... ignoring') # Check the UID in case this is a raw data file instead of an image file
		return
	seriesnum=im.SeriesNumber # Get the information to make the new file name
	seriesnum='%02d' % seriesnum
	seriesname=im.SeriesDescription
	seriesname=seriesname.replace('/','') # remove pesky forwardslashes!
	seriesname=seriesname.replace('.','') # remove pesky dots!
	seriesname=seriesname.replace('_','') # remove pesky underscores!
	seriesname=seriesname.replace(':','') # Remove pesky colons!
	imagenumber=im.InstanceNumber
	imagenumber='%04d' % imagenumber
	newname=str(seriesnum)+'_'+seriesname+'_'+str(imagenumber)+'.dcm' # add a .dcm extension so we can find them all again afterwards
	os.rename(filename,os.path.join(dstfolder,newname)) # Rename the file and move to destination

def dicom_renameSiemens(filename,dstfolder):

	import dicom
	import os

	im=dicom.read_file(filename) # Read the file
	seriesnum=im.SeriesNumber # Get the information to make the new file name
	seriesnum='%03d' % seriesnum
	if hasattr(im,'SeriesDescription'):
		seriesname=im.SeriesDescription
	elif hasattr(im, 'SequenceName'):
		seriesname=im.SequenceName
	else:
		seriesname='POSSIBLE_NON_IMAGE'

	seriesname=seriesname.replace('/','') # remove pesky forwardslashes!
	seriesname=seriesname.replace('.','') # remove pesky dots!
	seriesname=seriesname.replace('_','') # remove pesky underscores!
	seriesname=seriesname.replace(':','') # Remove pesky colons!
	seriesname=seriesname.replace('[','') # Remove pesky square brackets!
	seriesname=seriesname.replace(']','') # Remove pesky square brackets!
	imagenumber=im.InstanceNumber
	imagenumber='%04d' % imagenumber
	newname=str(seriesnum)+'_'+seriesname+'_'+str(imagenumber)+'.dcm' # add a .dcm extension so we can find them all again afterwards
	#newname=str(seriesnum)+'_'+str(imagenumber)+'_'+seriesname+'.ima' # add a .dcm extension so we can find them all again afterwards
	os.rename(filename,os.path.join(dstfolder,newname)) # Rename the file and move to destination

###############################################################################################
# Helper function to sort .dcm files into folders for each series after renaming as above

def sort_series(direct):
	import os
	import glob

	print('Sorting dicom files')
	allnames=glob.glob(os.path.join(direct,'*.dcm')) # find the dicom files
	allnames.sort()
	filenames=[os.path.split(X)[1] for X in allnames]# find the filenames only
	filenames.sort()
	allseries=[X.split('_')[0]+'_'+X.split('_')[1] for X in filenames]# find the series numbers and names
	# allseries=[X.split('_')[0] for X in filenames]
	series=set(allseries) # find unique entries in this list
	series=list(series) # convert back to a list for convenience

	# make directories with these names and move files with matching beginning of filename into them
	for Y in series:
		os.mkdir(os.path.join(direct,Y))
		#print(os.path.join(direct,(Y+'*')))
		#print(glob.glob(os.path.join(direct,(Y+'*.dcm'))))
		for Z in glob.glob(os.path.join(direct,(Y+'*.dcm'))):
			os.rename(Z,os.path.join(direct,Y,os.path.split(Z)[1]))

def sort_DWI():
	# Run this in the DWI series folder
	import os
	import glob
	import shutil

	filenames=glob.glob('*.dcm')
	filenames.sort()
	b800=filenames[2:240:3]
	b100=filenames[0:240:3]

	os.mkdir('b100')
	os.mkdir('b800')

	for x in b800:
		shutil.copy(x,'b800')
	for x in b100:
		shutil.copy(x,'b100')

	shutil.move('b100','..')
	shutil.move('b800','..')

###############################################################################################
# Helper function to rename folders for each series after renaming as above
def renamefolders(foldernames):	
	import os
	import glob
	for j in range(len(foldernames)):
		#Get the name of the first file:
		try:
			firstfile=glob.glob(os.path.join(foldernames[j],'*.dcm'))
			firstfile.sort()
			firstfile=firstfile[0]
			folderpath=os.path.split(foldernames[j])[0]
			newfoldername=os.path.join(folderpath,os.path.split(firstfile)[1].split('_0001')[0])
			os.rename(foldernames[j],newfoldername)
		except IndexError:
			pass

