#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob
sys.path.append('../nmutils')
from utils import *
import numpy as np
import matplotlib.pyplot as plt
import warnings
import h5py
warnings.filterwarnings("ignore")

### Input ###

PATH = '../resources/'
FILE = 'example.h5'
SAVE_PATH = '../resources/'

### Collecting data informations ###

try:
	inputFile = h5py.File(PATH + FILE,'r')
	print "File " + FILE + " loaded"
except ValueError:
	raise

### Creating matrix from data ###

rawData = np.array(inputFile['data/data'])
rawTheta = np.array(inputFile['data/theta'])

theta = np.sort(rawTheta)
sinogramData = np.zeros(np.shape(rawData))
reconstructedData = np.zeros((np.size(rawData,1),np.size(rawData,1),np.size(rawData,2)))

### Corrections ###

#Correcting unsorted 2theta acquisition
argsortVal = np.argsort(rawTheta)
if not np.array_equal(theta,rawTheta):
	sorting = 0
	while sorting < np.max(argsortVal)-1:
		sinogramData[sorting,:,:] = rawData[argsortVal[sorting],:,:]
		sorting = sorting+1
		progression("Sorting data................ ",np.size(argsortVal)-1,sorting)
print

#Removing outlier pixels from data
for i in range(0,np.size(rawData,2)):
	sinogramData[:,:,i] = findOutlierPixels(sinogramData[:,:,i],tolerance=3,worry_about_edges=True)	
	sinogramData[:,:,i] = divideByFirstColumn(sinogramData[:,:,i])
	progression("Correcting wrong pixels..... ",np.size(rawData,2),i)
print

#Correcting thermal/beam drifts
for i in range(0,np.size(rawData,2)):
	CoM = centerOfMass(sinogramData[:,:,i],axis=1)
	sinogramData[:,:,i] = fixDrift(sinogramData[:,:,i],CoM)
	progression("Correcting drifts........... ",np.size(rawData,2),i)
print

### Reconstruction ###

for i in range(0,np.size(rawData,2)):
	reconstructedData[:,:,i] = reconstruction(sinogramData[:,:,i],theta,output_size=np.size(rawData,1))
	progression("Reconstructing data......... ",np.size(rawData,2),i)
print

saveHdf5File(reconstructedData,SAVE_PATH,'reconstructedData.h5',mode='sliced')
saveImage(reconstructedData[:,:,25],SAVE_PATH,'test.png')
