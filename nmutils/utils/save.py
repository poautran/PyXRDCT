#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Importing couple of packages here
import numpy as np 
import sys, os, glob
import matplotlib.pyplot as plt
import h5py

#Definig all the functions required to save data

def saveImage(s,file_name):
	"""Saving image as a png image. File name as to be a string like 'image1'."""
	return scipy.misc.toimage(s).save('file_name'+'png')

def saveHdf5File(s,file_name,mode='stack'):
	"""Saving a matrix as a stack or as multiple images in a hdf5 file. Switch saving mode between 
	'stack' and 'sliced'. File name as to be a string like 'file1'"""
	f_out = h5py.File(file_name,'w')
	if mode=='stack':
		dset = im_fluo_out.create_dataset('data/'+file_name, np.shape(s), dtype='f')
		dset = s
	if mode=='sliced':
		for i in range(0,np.size(s,2)):
			dset = im_fluo_out.create_dataset('data/'+file_name+'%0d'%i, np.shape(s[:,:,0]), dtype='f')
			dset = s[:,:,i]
			sys.stdout.write("\r\x1b[K"+ "Saving data as hdf5 file: " + str( '%.2f'%float((float(i)+1)*100/np.size(s,2)))+ "%")	
			sys.stdout.flush()

