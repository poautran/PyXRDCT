#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob
PROGPATH = os.path.realpath(__file__)
PROGPATH = PROGPATH[:-17]
sys.path.append(PROGPATH + 'nmutils')
from utils import *
import numpy as np
import matplotlib.pyplot as plt
import warnings
import h5py
warnings.filterwarnings("ignore")
import argparse
import fabio
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import re

### Input ###

def main():
	parser = argparse.ArgumentParser(description='This program generated .xrdct files to be processed by the core. Input of the program must be EDF images (HDF will be supported soon). The program currently takes file names such as sample_STEP_ROT.edf (Ex: sample_0022_0011.edf (rotation number 11, step number 22))')
	parser.add_argument('INPUT',metavar="INPUT", help="List of files. Ex : *.edf", nargs='+')
	parser.add_argument('-j','-json',dest="JSON", default=".azimint.json", help="Configuration file containing the processing to be done")
	parser.add_argument('-o','--output',type=str,dest='OUTPUT',help='Output path')
	parser.add_argument('-m','--multiplier',type=int,dest='MULTIPLIER',help='Multiplies data by given number')
	parser.add_argument('-n','--normalize',help='Normalizes data from average at high angle (last 10%% -1 of the pattern)',dest='NORMALIZE',action='store_true')
	parser.add_argument('-t','--theta',help='Store a different value of 2 theta if acquisition is not in correct theta order (Ex -t 11,180 is the command to build theta as theta = i*11%%180)',dest='THETA')
	parser.set_defaults(func=run)
	args = parser.parse_args()
	args.func(args)

def run(args):
	FILE = args.INPUT
	if np.shape(FILE)==(1,):
		FILE = np.sort(glob.glob(str(FILE[0])))
	sample_name = str(raw_input("Enter sample name for saving: "))
	SAVE_PATH = os.getcwd()
	if args.OUTPUT:
		SAVE_PATH = args.OUTPUT
	else:
		print('!!! Warning files will be saved in the current folder because no output was defined.')
	currentFile = np.zeros((len(FILE),2))

	### Read json file ###
	jsonParam = readJson(args.JSON)

	### Grabbing first file to check matrix size ###
	for i in range(0,len(FILE)):
		currentFile[i] = re.findall(r'\d{3,7}',FILE[i])
		progression("Cheking files, matrix size definition ..... ",len(FILE),i)
	pattern = np.zeros((int(np.max(currentFile[:,0])),int(np.max(currentFile[:,1])+1),int(jsonParam['nbpt_rad'])))  

	### Integration of FILE ###
	azimutalIntegrator = AzimuthalIntegrator(dist=jsonParam['dist'], poni1=jsonParam['poni1'], poni2=jsonParam['poni2'], rot1=jsonParam['rot1'], rot2=jsonParam['rot2'], rot3=jsonParam['rot3'], pixel1=jsonParam['pixel1'], pixel2=jsonParam['pixel2'], splineFile=jsonParam['splineFile'], detector=jsonParam['detector'], wavelength=jsonParam['wavelength'])
	dark = np.array(fabio.open(jsonParam['dark_current']).data)
	flat = np.array(fabio.open(jsonParam['flat_field']).data)
	mask = np.array(fabio.open(jsonParam['mask_file']).data)
	offset_rot = int(np.min(currentFile[:,0]))
	offset_trans = int(np.min(currentFile[:,1]))
	for i in range(len(FILE)):
		dataFile = np.array(fabio.open(FILE[i]).data)
		dataX, dataY = azimutalIntegrator.integrate1d(dataFile, int(jsonParam['nbpt_rad']), filename=None, correctSolidAngle=jsonParam['do_solid_angle'], variance=None, error_model=None, radial_range=(float(jsonParam['radial_range_min']),float(jsonParam['radial_range_max'])), azimuth_range=None, mask=mask, dummy=jsonParam['do_dummy'], delta_dummy=jsonParam['delta_dummy'], polarization_factor=None, method='csr', dark=dark, flat=flat, unit=jsonParam['unit'], safe=True, normalization_factor=1, profile=False, all=False, metadata=None)
		#if args.SEPARATE:
		#	bragg, amorphous = azimutalIntegrator.separate(dataFile, npt_rad=1024, npt_azim=512, unit=jsonParam['unit'], method='splitpixel', percentile=50, mask=None, restore_mask=True)
		currentFile[i] = re.findall(r'\d{3,7}',FILE[i])
		pattern[int(currentFile[i,0])-offset_rot,int(currentFile[i,1])-offset_trans,:] = dataY
		progression("Integrating data............. ",len(FILE),i)
	print

	### Storing special 2-theta ###	
	theta = np.linspace(int(np.min(currentFile[:,0])),int(np.max(currentFile[:,0])),int(np.max(currentFile[:,0])))	
	if args.THETA:
		special_theta = np.fromstring(args.THETA,dtype=int,sep=',')
		for i in range(0,np.size(theta,0)):
			theta[i] = i*(int(special_theta[0])) % int(special_theta[1])

	### Normalize ###
	if args.NORMALIZE:
		normValue = int(np.around(np.size(pattern,2)*0.05, decimals=1))
		print normValue
		for i in range(0,np.size(pattern,0)):
			for j in range(0,np.size(pattern,1)):
				pattern[i,j,:] = pattern[i,j,:]/np.average(pattern[i,j,:-normValue])

	### Multiplier ###
	if args.MULTIPLIER:
		pattern = pattern*int(args.MULTIPLIER)

	### Saving ###
	if args.OUTPUT:
		f = h5py.File(args.OUTPUT + sample_name + '_sinogram.xrdct','w')
	else:
		f = h5py.File(sample_name + '_sinogram.xrdct','w')
	grp = f.create_group("data")

	dset = f.create_dataset('data/data', np.shape(pattern), dtype='f')
	dset[:,:,:] = pattern[:,:,:]

	dset_theta = f.create_dataset('data/theta', np.shape(theta), dtype='f')
	dset_theta[:] = theta[:]

	dset_dataX = f.create_dataset('data/dataX', np.shape(dataX), dtype='f')
	dset_dataX[:] = dataX[:]

if __name__=="__main__":
	main()

