#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Importing couple of packages here
import numpy as np 
import sys, os, glob
import matplotlib.pyplot as plt
import h5py
from display import *
from scipy import misc

#Definig all the functions required to save data

def saveImage(s,save_path,file_name):
	"""Saving image as a png image. File name as to be a string like 'image1'."""
	if not os.path.exists(save_path+'results/'):
		os.makedirs(save_path+'results/')
	misc.toimage(s).save(save_path+'results/'+file_name)
	print 'File '+file_name+ ' saved'

def saveHdf5File(s,save_path,file_name,mode='stack'):
	"""Saving a matrix as a stack or as multiple images in a hdf5 file. Switch saving mode between 
	'stack' and 'sliced'. File name as to be a string like 'file1'"""
	if not os.path.exists(save_path+'results/'):
		os.makedirs(save_path+'results/')
	f_out = h5py.File(save_path+'results/'+file_name,'w')
	if mode=='stack':
		dset = f_out.create_dataset('data/'+file_name[:-3], np.shape(s), dtype='f')
		dset[:,:,:] = s
	if mode=='sliced':
		for i in range(0,np.size(s,2)):
			dset = f_out.create_dataset('data/'+file_name[:-3]+'_%05d'%i, np.shape(s[:,:,0]), dtype='f')
			dset[:,:] = s[:,:,i]
			progression("Saving data as hdf5 file ",np.size(s,2),i)
		print
	print 'File '+file_name+ ' saved'