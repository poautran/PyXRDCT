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
import re

### Input ###

def main():
	parser = argparse.ArgumentParser(description='This program generated .xrdct files to be processed by the core. Input of the program must be ASCII two column files (.dat,.txt...). The program currently takes file names such as sample_STEP_ROT.dat (Ex: sample_0022_0011.dat (rotation number 11, step number 22))')
	parser.add_argument('INPUT',metavar="INPUT", help="List of files. Ex : *.dat", nargs='+')
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
	current_file = np.zeros((len(FILE),2))
	print 
	
	### Grabbing first file to check matrix size

	for i in range(0,len(FILE)):
		current_file[i] = re.findall(r'\d{3,7}',FILE[i])
		progression("Cheking files, matrix size definition ..... ",len(FILE),i)
	row_lines = np.genfromtxt(FILE[0],dtype=float,skip_header=23)
	pattern = np.zeros((int(np.max(current_file[:,0])),int(np.max(current_file[:,1])+1),int(np.size(row_lines,0))))

	### Starting sinogram stacking ###

	offset_rot = int(np.min(current_file[:,0]))
	offset_trans = int(np.min(current_file[:,1]))
	for i in range(0,len(FILE)):
		current_file[i] = re.findall(r'\d{3,7}',FILE[i])
		row_lines = np.genfromtxt(FILE[i],dtype=float,skip_header=23)
		pattern[int(current_file[i,0])-offset_rot,int(current_file[i,1])-offset_trans,:] = row_lines[:,1]
		progression("Stacking data................ ",len(FILE),i)
	
	### Storing special 2-theta ###	
	
	theta = np.linspace(int(np.min(current_file[:,0])),int(np.max(current_file[:,0])),int(np.max(current_file[:,0])))
	if args.THETA:
		special_theta = np.fromstring(args.THETA,dtype=int,sep=',')
		for i in range(0,np.size(theta,0)):
			theta[i] = i*(int(special_theta[0])) % int(special_theta[1])
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

if __name__=="__main__":
	main()

