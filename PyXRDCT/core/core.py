#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob
PROGPATH = os.path.realpath(__file__)
PROGPATH = PROGPATH[:-12]
sys.path.append(PROGPATH + 'nmutils')
from utils import *
import numpy as np
import matplotlib.pyplot as plt
import warnings
import h5py
warnings.filterwarnings("ignore")
import argparse
from scipy.ndimage import imread

### Input ###

def main():
	parser = argparse.ArgumentParser(description='This program correct and reconstruct XRD-CT data. See https://github.com/poautran/PyXRDCT for help.')
	parser.add_argument('INPUT',type=str,help='Input file (.xrdct from generate .xrdct (supported) or .h5 from pyFAI: diff_tomo. See https://pyfai.readthedocs.io/en/latest/man/diff_tomo.html) (not supported yet)')
	parser.add_argument('-o','--output',type=str,dest='OUTPUT',help='Output path')
	parser.add_argument('-d','--delete',help='Remove live in the sinogram from calculation. Usage: -d 180,179,114,103 to delete lines 180 179 114 and 103 ',dest='DELETE')
	parser.add_argument('-c','--CoM',help='Correct thermal drifts, beam drift using center of mass',dest='CORRECT',action='store_true')
	parser.add_argument('-a','--air',type=str,help='Corrects effect of air with input pattern',dest='AIR')
	parser.add_argument('-R','--overwrite',help='Overwrites the calculated sinogram with a new one',dest='OVERWRITE',action='store_true')
	parser.add_argument('-n','--normalize',help='Normalizes data from average at high angle',dest='NORMALIZE',action='store_true')
	parser.add_argument('-f','--filter',type=str,help='Multiplies the sinograms by a filter to remove contribution of the air. Input must be a binary image with 0 for air and 1 for sample',dest='FILTER')
	parser.add_argument('-ol','--outliers',help='Remove outliers of the sinogram',dest='OUTLIERS',action='store_true')
	parser.add_argument('-r','--reconstruct',help='Reconstruct data from corrected sinogram using Filtered Back Projection algorithm',dest='RECONSTRUCT',action='store_true')
	parser.add_argument('-p','--papyrus',type=str,help='Substract payrus signal ',dest='PAPYRUS')
	parser.set_defaults(func=run)
	args = parser.parse_args()
	args.func(args)

def run(args):
	FILE = args.INPUT
	FILE_NO_EXTENSION = FILE[:-6]
	SAVE_PATH = os.getcwd()
	print SAVE_PATH
	if args.OUTPUT:
		SAVE_PATH = args.OUTPUT
	else:
		print('!!! Warning files will be saved in the current folder because no output was defined.')

	REFERENCE_SLICE_NUMBER = 500

	### Collecting data informations ###

	try:
		inputFile = h5py.File(FILE,'r')
		print "File " + FILE + " loaded"
	except ValueError:
		raise

	### Creating matrix from data ###

	rawData = np.array(inputFile['data/data'])
	rawTheta = np.array(inputFile['data/theta'])

	theta = np.sort(rawTheta)
	rawData = rawData[:,:288,:]

	sinogramData = np.zeros(np.shape(rawData))
	reconstructedData = np.zeros((np.size(rawData,1),np.size(rawData,1),np.size(rawData,2)))

	### Corrections ###

	#Correcting unsorted 2theta acquisition
	argsortVal = np.argsort(rawTheta)
	if not np.array_equal(theta,rawTheta):
		sorting = 0
		while sorting < np.max(argsortVal)-1:
			sinogramData[sorting,:,:] = rawData[argsortVal[sorting],:,:]
			progression("Sorting data................ ",np.size(argsortVal)-2,sorting)
			sorting = sorting+1
	print

	#Normalizing data
	if args.NORMALIZE:
		for i in range(0,np.size(rawData,2)):
			sinogramData[:,:,i] = normalize(sinogramData[:,:,i])
			progression("Normalizing data............ ",np.size(rawData,2),i)
		print 

	##Deleting lines
	if args.DELETE:
		deleted_line = np.fromstring(args.DELETE,dtype=int,sep=',')
		for i in range(0,len(deleted_line)):
			sinogramData = np.delete(sinogramData, deleted_line[i], axis=0)
			theta = np.delete(theta, deleted_line[i], axis=0)
			progression("Deleting lines.............. ",len(deleted_line),i)
		print

	#Removing outlier pixels from data
	if args.OUTLIERS:
		for i in range(0,np.size(rawData,2)):
			sinogramData[:,:,i] = findOutlierPixels(sinogramData[:,:,i],tolerance=10,worry_about_edges=False)	
			progression("Correcting wrong pixels..... ",np.size(rawData,2),i)
		print

	#Correcting thermal/beam drifts
	if args.CORRECT:
		CoM = centerOfMass(sinogramData,axis=1,ref_slice=REFERENCE_SLICE_NUMBER)
		for i in range(0,np.size(rawData,2)):
			sinogramData[:,:,i] = fixDrift(sinogramData[:,:,i],CoM)
			progression("Correcting drifts........... ",np.size(rawData,2),i)
		print

	### Substract air from raw data ###
	if args.AIR:
		dataAir = np.genfromtxt(args.AIR,dtype=float)
		FILE_NO_EXTENSION = FILE_NO_EXTENSION + '_SUBAIR'
		for i in range(0,np.size(sinogramData,0)):
			for j in range(0,np.size(sinogramData,1)):
				currentAir = dataAir[:,1]*(sinogramData[i,j,96]/dataAir[96,1])
				sinogramData[i,j,:] = sinogramData[i,j,:]-currentAir
			progression("Substacting air............. ",np.size(sinogramData,0),i)
			plt.show()
		print
		

	### Filter ###
	if args.FILTER:
		filterData = args.FILTER
		filterImage = imread(filterData)
		for i in range(0,len(deleted_line)-1):
			filterImage = np.delete(filterImage, deleted_line[i+1], axis=0)
		for i in range(0,np.size(rawData,2)):
			sinogramData[:,:,i] = sinogramData[:,:,i]*filterImage
			progression("Masking air................. ",np.size(rawData,2),i)
		print

	### Saving ###
	if (args.OVERWRITE == True or os.path.isfile(FILE_NO_EXTENSION+'_corrected.h5') == False):
		saveHdf5File(sinogramData,SAVE_PATH,FILE_NO_EXTENSION+'_corrected.h5',mode='sliced')
	else:
		print('!!! Warning sinogram file exists, use command -R to overwrite it')

	### Reconstruction ###
	if args.RECONSTRUCT:
		for i in range(0,np.size(rawData,2)):
			reconstructedData[:,:,i] = reconstruction(sinogramData[:,:,i],theta,output_size=np.size(rawData,1))
			progression("Reconstructing data......... ",np.size(rawData,2),i)
		print
		
	
	### Substract Papyrus ###
	if args.PAPYRUS:
		dataPapyrus = np.genfromtxt(args.PAPYRUS,dtype=float)
		FILE_NO_EXTENSION = FILE_NO_EXTENSION + '_SUBPAP'
		for i in range(0,np.size(reconstructedData,0)):
			for j in range(0,np.size(reconstructedData,1)):
				if np.max(reconstructedData[i,j,:]) >= 0.45:
					reconstructedData[i,j,:] = reconstructedData[i,j,:] - (dataPapyrus[:,1]*(reconstructedData[i,j,REFERENCE_SLICE_NUMBER]/dataPapyrus[REFERENCE_SLICE_NUMBER,1]))

		
	if args.OVERWRITE:
		saveHdf5File(reconstructedData,SAVE_PATH,FILE_NO_EXTENSION+'_reconstructed_stack.h5',mode='stack')
	else:
		print('!!! Warning reconstruction file exists, use command -R to overwrite it')


if __name__=="__main__":
	main()



