# Helper funditon to write GSK spreadsheet for sjogrens data

import csv

# write to csv, each result on a new line
def send_to_file(subjectID,res_line, Map, region, SEQ,UMCU=0):
	#res_line is [mean, sd, median, lq, uq]
	#SEQ is the sequence variable that has to be incremented for each result written out
	# Set up variables to write results
	STUDYID=203818
	SUBJID=subjectID
	VISITNUM=30
	VISIT='VISIT 1'
	PTMNUM=140
	PTM='VISIT 1'
	TPTREFN=70
	TPTREF='VISIT 1'
	ELTMNUM=''
	ELTMUNIT=''
	ACTTM=''
	EXMNUM=''
	MRIRESTX=''
	MRITYPECD=12

	###### Make dictionaries
	#Region labels
	REGINTCDdict={'lt_parotid':536,'rt_parotid':537,'lt_submandibular':538,'rt_submandibular':539,'parotid':547, 'submandibular':548}

	#Scan dates and times labels
	SCNSTDTdict={
			'101':['14-Dec-16','14-Dec-16','13:32','14:15'],
			'102':['08-Dec-16','08-Dec-16','11:06','11:49'],
			'104':['21-Aug-17','21-Aug-17','17:14','17:58'],
			'106':['07-Aug-17','07-Aug-17','13:08','13:48'],
			'108':['26-Feb-18','26-Feb-18','10:20','11:02'],
			'111':['05-Mar-18','05-Mar-18','10:11','10:55'],
			'113':['05-Mar-18','05-Mar-18','12:20','13:06'],
			'114':['23-Mar-18','23-Mar-18','11:28','12:10'],
			'115':['21-Mar-18','21-Mar-18','10:42','11:25'],
			'116':['05-Apr-18','05-Apr-18','11:55','12:44'],
			'117':['30-Apr-18','30-Apr-18','11:28','12:11'],
			'118':['10-May-18','10-May-18','10:43','11:31'],
			'401':['16-Mar-17','16-Mar-17','14:18','14:44'],
			'402':['27-Mar-17','27-Mar-17','14:45','15:32'],
			'403':['30-May-17','30-May-17','13:14','14:07'],
			'404':['12-Jul-17','12-Jul-17','13:11','13:56'],
			'405':['20-Oct-17','20-Oct-17','10:13','11:10'],
			'406':['06-Feb-18','06-Feb-18','09:51','10:28'],
			'407':['23-Apr-18','23-Apr-18','14:10','14:52'],
			'408':['18-Apr-18','18-Apr-18','11:26','12:13'],
			'600':['18-Apr-17','18-Apr-17','12:34','13:15'],
			'601':['12-Jun-17','12-Jun-17','13:20','14:06'],
			'602':['28-Oct-17','28-Oct-17','10:27','11:11'],
			'603':['30-Oct-17','30-Oct-17','11:01','11:38']
		}

	#Value code labels - order is [mean, sd, median, lq, uq, UNITS]
	MRITSTCDdict={
			'T1map':[155,141,140,142,143,'MSEC'],
			'MaxEnhancement':[154,137,136,138,139,'%'],
			'MaxEnhancementConc':[154,137,136,138,139,'MMOL'],
			'InitialRateEnhancement':[153,133,132,134,135,'%/SEC'],
			'InitialRateEnhancementConc':[153,133,132,134,135,'MMOL/SEC'],
			'ExtKetyfitConc':[152,145,144,146,147,'/MIN'],
			'ExtKetyfitSI':[152,145,144,146,147,'/MIN'],
		}
	#######
	#Choose dictionary values for this result

	REGINTCD=REGINTCDdict[region]
	SCNSTDT=SCNSTDTdict[SUBJID]
	MRITSTCDarray=MRITSTCDdict[Map]
	MRIRESU=MRITSTCDarray[-1]
	if UMCU==1:
		SUBJID=SUBJID[0:4]# Now remove 'followup' if UMCU data

	valuearray=['mean','sd','median','lq','uq']

	with open ('/Users/lkershaw/Desktop/newresults.csv','a',newline='') as csvfile:
		resultswriter=csv.writer(csvfile,delimiter=',')

		for i in range(5):
			value=valuearray[i]
			MRITSTCD=MRITSTCDarray[i]
			MRIRESN=res_line[i]
			resultswriter.writerow([STUDYID,'100'+SUBJID,VISITNUM,VISIT,PTMNUM,PTM,TPTREFN,TPTREF,ELTMNUM,ELTMUNIT,SEQ+i,ACTTM,EXMNUM,SCNSTDT[0],SCNSTDT[1],SCNSTDT[2],SCNSTDT[3],MRITSTCD,MRIRESN,MRIRESTX,REGINTCD,MRITYPECD,MRIRESU])


