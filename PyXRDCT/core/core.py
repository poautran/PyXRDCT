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
from diffpy.pdfgetx import PDFConfig
from diffpy.pdfgetx import PDFGetter

### Input ###

def main():
	parser = argparse.ArgumentParser(description='This program correct and reconstruct XRD-CT data. See https://github.com/poautran/PyXRDCT for help.')
	parser.add_argument('INPUT',type=str,help='Input file (.xrdct from generate .xrdct (supported) or .h5 from pyFAI: diff_tomo. See https://pyfai.readthedocs.io/en/latest/man/diff_tomo.html) (not supported yet)')
	parser.add_argument('-o','--output',type=str,dest='OUTPUT',help='Output path')
	parser.add_argument('-d','--delete',help='Remove live in the sinogram from calculation. Usage: -d 180,179,114,103 to delete lines 180 179 114 and 103 ',dest='DELETE')
	parser.add_argument('-c','--CoM',help='Correct thermal drifts, beam drift using center of mass',dest='CORRECT',action='store_true')
	parser.add_argument('-a','--air',type=str,help='Corrects effect of air with input pattern',dest='AIR')
	parser.add_argument('-e','--extra',type=str,help='Corrects effect of extra contribution with input pattern',dest='EXTRA')
	parser.add_argument('-p','--pdf',type=str,help='Extract PDF signal from sinogram',dest='PDF')
	parser.add_argument('-R','--overwrite',help='Overwrites the calculated sinogram with a new one',dest='OVERWRITE',action='store_true')
	parser.add_argument('-ol','--outliers',help='Remove outliers of the sinogram',dest='OUTLIERS',action='store_true')
	parser.add_argument('-r','--reconstruct',help='Reconstruct data from corrected sinogram using Filtered Back Projection algorithm',dest='RECONSTRUCT',action='store_true')
	parser.set_defaults(func=run)
	args = parser.parse_args()
	args.func(args)

def run(args):
	FILE = args.INPUT
	FILE_NO_EXTENSION = FILE[:-6]
	SAVE_PATH = os.getcwd()
	print(SAVE_PATH)
	if args.OUTPUT:
		SAVE_PATH = args.OUTPUT
	else:
		print('!!! Warning files will be saved in the current folder because no output was defined.')

	REFERENCE_SLICE_NUMBER = 200

	### Collecting data informations ###

	try:
		inputFile = h5py.File(FILE,'r')
		print("File " + FILE + " loaded")
	except ValueError:
		raise

	### Creating matrix from data ###

	rawData = np.array(inputFile['data/data'])
	rawTheta = np.array(inputFile['data/theta'])
	dataX = np.ndarray.flatten(np.array(inputFile['data/dataX']))

	theta = np.sort(rawTheta)
	#rawData = rawData[:,:288,:]

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
	print()

	##Deleting lines
	if args.DELETE:
		deleted_line = np.fromstring(args.DELETE,dtype=int,sep=',')
		for i in range(0,len(deleted_line)):
			sinogramData = np.delete(sinogramData, deleted_line[i], axis=0)
			theta = np.delete(theta, deleted_line[i], axis=0)
			progression("Deleting lines.............. ",len(deleted_line),i)
		print()

	### Removing outlier pixels from data ###
	if args.OUTLIERS:
		for i in range(0,np.size(rawData,2)):
			sinogramData[:,:,i] = findOutlierPixels(sinogramData[:,:,i],tolerance=10,worry_about_edges=False)	
			progression("Correcting wrong pixels..... ",np.size(rawData,2),i)
		print()

	### Subtract air from raw data ###
	if args.AIR:
		dataAir = np.genfromtxt(args.AIR,dtype=float)
		FILE_NO_EXTENSION = FILE_NO_EXTENSION + '_SUBAIR'
		for i in range(0,np.size(sinogramData,0)):
			for j in range(0,np.size(sinogramData,1)):
				currentAir = dataAir[:,1]*(sinogramData[i,j,96]/dataAir[96,1])
				sinogramData[i,j,:] = sinogramData[i,j,:]-currentAir
			progression("Substacting air............. ",np.size(sinogramData,0),i)
			plt.show()
		print()


	### Correcting thermal/beam drifts ###
	if args.CORRECT:
		CoM = centerOfMass(sinogramData,axis=1,ref_slice=REFERENCE_SLICE_NUMBER)
		for i in range(0,np.size(rawData,2)):
			sinogramData[:,:,i] = fixDrift(sinogramData[:,:,i],CoM)
			progression("Correcting drifts........... ",np.size(rawData,2),i)
		print()

	### Subtract extra pattern from raw data ###
	if args.EXTRA:
		dataExtra = np.genfromtxt(args.EXTRA,dtype=float)
		FILE_NO_EXTENSION = FILE_NO_EXTENSION + '_SUBPAP'
		for i in range(0,np.size(sinogramData,0)):
			for j in range(0,np.size(sinogramData,1)):
				currentExtra = dataExtra[:,1]*(sinogramData[i,j,147]/dataExtra[147,1])
				sinogramData[i,j,:] = sinogramData[i,j,:]-currentExtra
			progression("Substacting extra........... ",np.size(sinogramData,0),i)
		print()

	### Extract PDF signal ###
	sinogramDataPdf = np.copy(sinogramData)
	if args.PDF:
		cfg = PDFConfig()
		cfg.readConfig(args.PDF)
		pdfget = PDFGetter()
		pdfget.configure(cfg)
		sinogramDataPdf = np.zeros((np.size(sinogramData,0),np.size(sinogramData,0),round((cfg.rmax-cfg.rmin)/cfg.rstep)+1))
		for i in range(0,np.size(sinogramData,0)):
			for j in range(0,np.size(sinogramData,1)):
				currentPdfDataY = sinogramData[i,j,:]
				pdfget.getTransformation('gr')
				pdfget(dataX,currentPdfDataY)
				pdfResults = pdfget.results
				pdfResults = pdfResults[8]
				sinogramDataPdf[i,j,:] = pdfResults[1]
			progression("Extracting PDF.............. ",np.size(sinogramData,0),i)
		sinogramData = np.copy(sinogramDataPdf)
		FILE_NO_EXTENSION = FILE_NO_EXTENSION + '_PDF'
		print()

	### Saving ###
	if (args.OVERWRITE == True or os.path.isfile(FILE_NO_EXTENSION+'_corrected.h5') == False):
		saveHdf5File(sinogramData,SAVE_PATH,FILE_NO_EXTENSION+'_corrected.h5',mode='sliced')
	else:
		print('!!! Warning sinogram file exists, use command -R to overwrite it')

	### Reconstruction ###
	if args.RECONSTRUCT:
		for i in range(0,np.size(sinogramData,2)):
			reconstructedData[:,:,i] = reconstruction(sinogramData[:,:,i],theta,output_size=np.size(rawData,1))
			progression("Reconstructing data......... ",np.size(rawData,2),i)
		print()
		
	if args.OVERWRITE:
		saveHdf5File(reconstructedData,SAVE_PATH,FILE_NO_EXTENSION+'_reconstructed_stack.h5',mode='stack')
	else:
		print('!!! Warning reconstruction file exists, use command -R to overwrite it')


if __name__=="__main__":
	main()


