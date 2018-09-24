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
REFERENCE_SLICE_NUMBER = 20

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
<<<<<<< HEAD
		progression("Sorting data................ ",np.size(argsortVal)-1,sorting)
=======
>>>>>>> 5645e75bdaac69f1ed8885878ab503152bb61cc2
		sorting = sorting+1
		progression("Sorting data................ ",np.size(argsortVal)-1,sorting)
print

#Normalizing data
for i in range(0,np.size(rawData,2)):
	sinogramData[:,:,i] = normalize(sinogramData[:,:,i])
	sinogramData[:,:,i] = divideByFirstColumn(sinogramData[:,:,i])	
	progression("Normalizing data............ ",np.size(rawData,2),i)
print 

##Deleting lines
#deleted_line = []
#for i in range(0,np.size(deleted_line,0)-1):
#	sinogramData = np.delete(sinogramData, deleted_line[i], axis=0)
#	theta = np.delete(theta, deleted_line[i], axis=0)
#	progression("Deleting lines.............. ",np.size(deleted_line,0)-1,i)
#print	

#Correcting thermal/beam drifts
CoM = centerOfMass(sinogramData,axis=1,slice=REFERENCE_SLICE_NUMBER)
for i in range(0,np.size(rawData,2)):
	sinogramData[:,:,i] = fixDrift(sinogramData[:,:,i],CoM)
	progression("Correcting drifts........... ",np.size(rawData,2),i)
print

#Removing outlier pixels from data
for i in range(0,np.size(rawData,2)):
	sinogramData[:,:,i] = findOutlierPixels(sinogramData[:,:,i],tolerance=2,worry_about_edges=False)	
	progression("Correcting wrong pixels..... ",np.size(rawData,2),i)
print 

### Reconstruction ###

for i in range(0,np.size(rawData,2)):
	reconstructedData[:,:,i] = reconstruction(sinogramData[:,:,i],theta,output_size=np.size(rawData,1))
	progression("Reconstructing data......... ",np.size(rawData,2),i)
print

saveHdf5File(reconstructedData,SAVE_PATH,'reconstructedData.h5',mode='sliced')
saveImage(reconstructedData[:,:,25],SAVE_PATH,'test.png')





