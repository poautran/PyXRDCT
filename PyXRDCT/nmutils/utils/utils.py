#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Importing couple of packages here
import numpy as np 
import sys, os, glob
import matplotlib.pyplot as plt
import h5py
from skimage.transform import iradon 
from scipy.ndimage import shift
from scipy.ndimage.measurements import center_of_mass

#Definig all the functions required to correct and reconstruct sinograms

def findOutlierPixels(data,tolerance=3,worry_about_edges=True):
    """This function finds the hot or dead pixels in a 2D dataset. 
    tolerance is the number of standard deviations used to cutoff the hot pixels
    If you want to ignore the edges and greatly speed up the code, then set
    worry_about_edges to False.
    The function returns a list of hot pixels and also an image with with hot pixels removed"""

    from scipy.ndimage import median_filter
    blurred = median_filter(data, size=2)
    difference = data - blurred
    threshold = 10*np.std(difference)

    #Find the hot pixels, but ignore the edges
    hot_pixels = np.nonzero((np.abs(difference[1:-1,1:-1])>threshold) )
    hot_pixels = np.array(hot_pixels) + 1 #because we ignored the first row and first column

    fixed_image = np.copy(data) #This is the image with the hot pixels removed
    for y,x in zip(hot_pixels[0],hot_pixels[1]):
        fixed_image[y,x]=blurred[y,x]

    if worry_about_edges == True:
        height,width = np.shape(data)

        ###Now get the pixels on the edges (but not the corners)###

        #left and right sides
        for index in range(1,height-1):
            #left side:
            med  = np.median(data[index-1:index+2,0:2])
            diff = np.abs(data[index,0] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[index],[0]]  ))
                fixed_image[index,0] = med

            #right side:
            med  = np.median(data[index-1:index+2,-2:])
            diff = np.abs(data[index,-1] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[index],[width-1]]  ))
                fixed_image[index,-1] = med

        #Then the top and bottom
        for index in range(1,width-1):
            #bottom:
            med  = np.median(data[0:2,index-1:index+2])
            diff = np.abs(data[0,index] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[0],[index]]  ))
                fixed_image[0,index] = med

            #top:
            med  = np.median(data[-2:,index-1:index+2])
            diff = np.abs(data[-1,index] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[height-1],[index]]  ))
                fixed_image[-1,index] = med

        ###Then the corners###

        #bottom left
        med  = np.median(data[0:2,0:2])
        diff = np.abs(data[0,0] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[0],[0]]  ))
            fixed_image[0,0] = med

        #bottom right
        med  = np.median(data[0:2,-2:])
        diff = np.abs(data[0,-1] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[0],[width-1]]  ))
            fixed_image[0,-1] = med

        #top left
        med  = np.median(data[-2:,0:2])
        diff = np.abs(data[-1,0] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[height-1],[0]]  ))
            fixed_image[-1,0] = med

        #top right
        med  = np.median(data[-2:,-2:])
        diff = np.abs(data[-1,-1] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[height-1],[width-1]]  ))
            fixed_image[-1,-1] = med
    return fixed_image

def normalize(matrix):
    matrix_max, matrix_min = matrix.max(), matrix.min()
    matrix_normalized = (matrix - matrix_min)/(matrix_max - matrix_min)
    return matrix_normalized

def divideByFirstColumn(matrix):
    """This function devide a matrix by its first column to resolve
    wrong intemsity problems"""
    result = (matrix.T / matrix.sum(axis=1)).T
    return result

def centerOfMass(matrix,axis=1,ref_slice=0):
    """Calculating center of Mass of a matrix along column axis"""
    matrix_flat = matrix[:,:,ref_slice]
    x = np.linspace(0,np.size(matrix_flat,axis)-1,np.size(matrix_flat,axis))
    CoM = np.linspace(0,np.size(matrix_flat,abs(axis-1))-1,np.size(matrix_flat,abs(axis-1)))
    for i in range(0,np.size(matrix_flat,abs(axis-1))):
        CoM[i] = np.sum( (x*matrix_flat[i,:])/np.sum(matrix_flat[i,:]) )
    CoM = CoM - np.size(matrix_flat,axis)/2
    return CoM

def fixDrift(s,CoM):
    """Fixing beam, thermal drifts during acquisition. Assumes background subtracted also normalises."""
    sOut = np.zeros(np.shape(s))
    for i in range(0,np.size(s,0)):
        corrDrift = CoM[i]
        shift(s[i,:],-corrDrift,sOut[i,:],mode='constant')
    return sOut

def reconstruction(sinogram,theta,method='FBP',output_size='None',filter='ramp'):
	"""Reconstructing with Filtered Back Projection by deflaut. Needs theta as 1D matrix, output_size 
	as a number or default and filter could be ramp, cosine..."""
	if method=='FBP':
		if output_size=='default':
			rec = iradon(sinogram.T,theta,filter=filter)
		else:
			rec = iradon(sinogram.T,theta,output_size=output_size,filter=filter)
	return rec
