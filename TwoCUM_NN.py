# Module to make a neural network for fitting the 2CUM model

import torch
import torch.nn as nn
import torchvision.transforms as transforms
import torch.optim
import torchvision.datasets
import torch.nn.functional as F 
import numpy as np
import matplotlib.pyplot as plt
import os
import TwoCUM
import itertools

class training_data(object): # Training data class, which will be populated with... 
	def __init__(self):
		self.curves=[] #Simulated curves,
		self.AIF=[] #generated from this (these) AIF(s)
		self.AIF_source=[] #which are from here
		self.params=[] #using this set of parameters

	def read_AIF(self, AIF_path, tres): # Also input the time resolution for this AIF here (in s)
		self.AIF_source=self.AIF_source+list([AIF_path]) #Add the new path to the list of sources
		
		print('looking for AIF file...') #Read the AIF at that location
		if not os.path.isfile(AIF_path):
			print('AIF file not found')
			return	
		print('Found')
		AIFnew=np.squeeze(np.load(os.path.join(AIF_path)))
		if sum(self.AIF)>0:
			self.AIF=np.concatenate((self.AIF,AIFnew[:,np.newaxis]),1)
		else:
			self.AIF=AIFnew[:,np.newaxis]

		self.t=np.arange(0,tres*len(AIFnew),tres) #Make time points array

	def set_params(self, Enum, Fprange, vpnum, toffnum): #Set the possible ranges for the model parameters. Fprange=(mean, sd, size)
		#Aim here is to form a matrix of rows of E, Fp, vp, toff values.  Uniform distribution for E, vp, toff, normal for Fp

		Es=np.random.uniform(0,1,Enum)
		Fps=np.random.normal(Fprange[0],Fprange[1],Fprange[2])
		vps=np.random.uniform(0,1,vpnum)
		toffs=np.random.choice(self.t[0:10],toffnum) #random sample from first 10 time points

		self.params=np.array(list(itertools.product(Es,Fps,vps,toffs)))

	def calc_curve(self,params_index,AIF_index): #method to calculate a concentration curve from a given parameter set index and AIF index
		conc_curve=TwoCUM.TwoCUM(self.params[params_index,0:3],self.t,np.squeeze(self.AIF[:,AIF_index]),self.params[params_index,3])
		self.conc_curve=conc_curve
		return conc_curve

class Net(nn.Module):
	def __init__(self):
		super(Net, self).__init__()
		self.fc1=nn.Linear(150,150) #150 in input, 150 in first hidden
		self.fc2=nn.Linear(150,75) #150 from first hidden, 75 in second hidden
		self.fc3=nn.Linear(75,4) #75 from second hidden, 4 output
	
	def forward(self,x):
		x=F.relu(self.fc1(x))
		#print(x)
		x=F.relu(self.fc2(x))
		#print(x)
		x=F.relu(self.fc3(x))
		#print(x)	
		return x






