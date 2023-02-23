#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Project: PyXRDCT
#             https://github.com/poautran/PyXRDCT
#
#    Copyright (C) 2022-2023 European Synchrotron Radiation Facility, Grenoble,
#             France
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
import h5py, hdf5plugin
import os, sys, time

import readh5
import slurm
import saveh5

from skimage.transform import iradon

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('image', cmap='gray')

def shift_sino(s,s_shift):
    from scipy.ndimage import shift
    sOut = np.zeros(np.shape(s))
    for i in range(0,np.size(s,1)):
        shift(s[:,i],-s_shift,sOut[:,i],mode='nearest')
    return sOut

class Reconstruction:
    """
    Initialise reconstruction class
    """
    def __init__(self, readH5Input, intFile=False):
        self.data = readH5Input
        if intFile:
            self.integrate = intFile

    def reconstruct2d_s3dxrd(self,  binning=1, shift=0, plot=False, save=True):
        """
        Reconstructs 2D slice of grains from segmented s3DXRD.
        """
        tdxrdData = np.empty((len(self.data.y),len(self.data.rot[0])), dtype=np.float32)
        with h5py.File(os.path.join(self.data.savePath,self.data.dataset+'_s3dxrd_segmented.h5'),'r') as h5In:
            for i,scan in enumerate(self.data.scans): 
                tdxrdData[i] = h5In[scan]['nnz'][:]
        tdxrdDataSino, a,y = np.histogram2d( np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),weights = np.array(tdxrdData).ravel(), bins=(int(self.data.rot.shape[1]/binning),int(self.data.rot.shape[0]/binning)) )
        tdxrdDataRecon = iradon(shift_sino(tdxrdDataSino.T,shift),self.data.rot[0],circle = True ,output_size = int(len(self.data.y)/binning))
        if save:
            saveh5.saveReconstructedH5(os.path.join(self.data.savePath,self.data.dataset+'_s3dxrd_2dreconstruction.h5'),tdxrdDataRecon)
        if plot:
            plt.figure(figsize=(20,10))
            plt.subplot(121)
            plt.imshow(tdxrdDataSino)
            plt.title('%s: Segmented grains sinogram'%self.data.dataset)
            plt.subplot(122)
            plt.imshow(tdxrdDataRecon)
            plt.title('%s: Segmented grains reconstruction'%self.data.dataset)
            
    def reconstruct2d_xrdct(self, tths=[3,4], width=0.05, binning=1, shift=0, plot=False, save=True):
        """
        Reconstructs 2D slice of XRD-CT from provided array of energies.
        """
        TTH_MIN= 0.22
        TTH_MAX= 22
        xrdDataReconSave = []
        for tth in tths:
            idx = (np.abs(np.linspace(TTH_MIN,TTH_MAX,self.data.nbpt_rad) - tth)).argmin()
            idxWidth = int(((self.data.channels/(TTH_MAX-TTH_MIN))*width)
            xrdData = np.empty((len(self.data.y),len(self.data.rot[0])), dtype=np.float32)
            #Add integrated data
            with h5py.File(self.data.dataPath,'r') as h5In:
                for i,scan in enumerate(self.data.scans):
                    xrdData[i] = np.average(h5In[scan]['measurement'][self.data.xrddetector][:,idx-idxWidth:idx+idxWidth],axis=1)
            
                           
            xrdDataSino, a,y = np.histogram2d( np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),weights = np.array(xrdData).ravel(), bins=(int(self.data.rot.shape[1]/binning),int(self.data.rot.shape[0]/binning)) )
            xrdDataRecon = iradon(shift_sino(xrdDataSino.T,shift),self.data.rot[0],circle = True ,output_size = int(len(self.data.y)/binning))
            if plot:
                plt.figure(figsize=(20,10))
                plt.subplot(121)
                plt.imshow(xrdDataSino)
                plt.title('%s: XRD sinogram %skeV +/-%skeV'%(self.data.dataset,energy,width))
                plt.subplot(122)
                plt.imshow(xrdDataRecon)
                plt.title('%s: XRD reconstruction %skeV +/-%skeV'%(self.data.dataset,energy,width))
            xrdDataReconSave.append(xrdDataRecon)
        xrdDataReconSave = np.array(xrdDataReconSave)
        if save:
            saveh5.saveReconstructedH5(os.path.join(self.data.savePath,self.data.dataset+'_xrd_2dreconstruction.h5'),xrdDataReconSave,tths)
        
    def reconstruct2d_xrfct(self, energies=[2.013,28.612,49.127,61.140,0.5249,0.0543,4.952,28.612,49.127], binning=1, width=0.05, shift=0, plot=False, save=True):
        """
        Reconstructs 2D slice of XRF-CT from provided array of energies.
        """
        ENERGY_MAX = 81.92
        xrfDataReconSave = []
        for energy in energies:
            idx = (np.abs(np.linspace(0,ENERGY_MAX,self.data.channels) - energy)).argmin()
            idxWidth = int((self.data.channels/ENERGY_MAX)*width)
            xrfData = np.empty((len(self.data.y),len(self.data.rot[0])), dtype=np.float32)
            with h5py.File(self.data.dataPath,'r') as h5In:
                for i,scan in enumerate(self.data.scans):
                    xrfData[i] = np.average(h5In[scan]['measurement'][self.data.xrfdetector][:,idx-idxWidth:idx+idxWidth],axis=1)
            xrfDataSino, a,y = np.histogram2d( np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),weights = np.array(xrfData).ravel(), bins=(int(self.data.rot.shape[1]/binning),int(self.data.rot.shape[0]/binning)) )
            xrfDataRecon = iradon(shift_sino(xrfDataSino.T,shift),self.data.rot[0],circle = True ,output_size = int(len(self.data.y)/binning))
            if plot:
                plt.figure(figsize=(20,10))
                plt.subplot(121)
                plt.imshow(xrfDataSino)
                plt.title('%s: XRF sinogram %skeV +/-%skeV'%(self.data.dataset,energy,width))
                plt.subplot(122)
                plt.imshow(xrfDataRecon)
                plt.title('%s: XRF reconstruction %skeV +/-%skeV'%(self.data.dataset,energy,width))
            xrfDataReconSave.append(xrfDataRecon)
        xrfDataReconSave = np.array(xrfDataReconSave)
        if save:
            saveh5.saveReconstructedH5(os.path.join(self.data.savePath,self.data.dataset+'_xrf_2dreconstruction.h5'),xrfDataReconSave,energies)
            