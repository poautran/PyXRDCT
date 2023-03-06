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

import PyXRDCT.nmutils.utils.saveh5 as saveh5

class Input:
    """
    Class to input data for processing XRD/XRF/PDF-CT and 3DXRD.
    """
    def __init__(self, dataPath, session='Default'):
        """
        Initialize data path.
        """
        self.dataPath = dataPath
        self.dataset = os.path.basename(os.path.dirname(self.dataPath))
        self.sample = os.path.basename(os.path.dirname(os.path.dirname(self.dataPath)))
        self.expPath = os.path.dirname(os.path.dirname(os.path.dirname(self.dataPath)))
        if session='Default':
            self.savePath = os.path.join(self.expPath,'xrdct_analysis',self.sample,self.dataset)
        else:
            self.savePath = os.path.join(os.path.dirname(self.expPath),session,'xrdct_analysis',self.sample,self.dataset)
        saveh5.makeSaveDirs(self.savePath)

    def getScanGeometry(self):
        """
        Reads scan Geometry from h5 file.
        """
        scanKey=['1']
        scansInput = []
        self.y = []
        self.rot = []
        self.scans = []
        yMotor = self.yMotor()
        rotMotor = self.rotMotor()
        with h5py.File(self.dataPath,'r') as h5In:
            scans = list(h5In.keys())
            for scan in scans:
                if scan.split('.')[1] == scanKey[0]:
                    scansInput.append(int(scan.split('.')[0]))
        scansInput = sorted(scansInput)
        for scan in scansInput:
            self.scans.append('%s.%s'%(scan,scanKey[0]))
        with h5py.File(self.dataPath,'r') as h5In:
            for scan in self.scans:
                self.y.append(h5In['%s/instrument/positioners/%s'%(scan,yMotor)][()])
                self.rot.append(h5In['%s/measurement/%s'%(scan,rotMotor)][()])
        self.y = np.array(self.y)
        self.rot = np.array(self.rot)
        if len(np.array(self.y).shape) == 1:
            bufferY = np.copy(self.rot)
            for i in range(np.array(self.rot).shape[1]):
                bufferY[:,i] = self.y
            self.y = bufferY
        print('[INFO] Read %s scans'%max(scansInput))

    def yMotor(self):
        """
        Finds y motor from provided list.
        """
        yMotors = ['dty','diffty','diffy']
        with h5py.File(self.dataPath,'r') as h5In:
            for motor in yMotors:
                if motor in list(h5In['1.1/instrument/positioners'].keys()):
                    self.yMotor = motor
        print('[INFO] Y motor detected: %s'%self.yMotor)
        return self.yMotor

    def rotMotor(self):
        """
        Finds y motor from provided list.
        """
        rotMotors = ['rot','diffrz','shrz']
        with h5py.File(self.dataPath,'r') as h5In:
            for motor in rotMotors:
                if motor in list(h5In['1.1/measurement'].keys()):
                    self.rotMotor = motor
        print('[INFO] Rot motor detected: %s'%self.rotMotor)
        return self.rotMotor

    def getXrdDetector(self):
        """
        Finds XRD detector.
        """
        xrdDetectors = ['eiger','frelon3']
        with h5py.File(self.dataPath,'r') as h5In:
            for detector in xrdDetectors:
                if detector in list(h5In['1.1/measurement'].keys()):
                    self.xrddetector = detector
                    print('[INFO] XRD detector: %s'%self.xrddetector)

    def getXrfDetector(self):
        """
        Finds XRF detector and channels.
        """
        xrfDetectors = ['mca_det0']
        with h5py.File(self.dataPath,'r') as h5In:
            for detector in xrfDetectors:
                if detector in list(h5In['1.1/measurement'].keys()):
                    self.xrfdetector = detector
                    print('[INFO] XRF detector: %s'%self.xrfdetector)
                    self.channels = h5In['1.1/measurement'][detector][:].shape[1]
                    print('[INFO] XRF detector channels detected: %s'%self.channels)

    def getBeamMonitor(self):
        """
        Finds XRF detector and channels.
        """
        monitors = ['fpico6','fpico4']
        with h5py.File(self.dataPath,'r') as h5In:
            for monitor in monitors:
                if monitor in list(h5In['1.1/measurement'].keys()):
                    self.beamMonitor = monitor
                    print('[INFO] Beam Monitor: %s'%self.beamMonitor)

    def loadData(self):
        """
        Loads data Urls from bliss HDF5.
        """
        self.getXrfDetector()
        self.getXrdDetector()
        self.getBeamMonitor()
        self.getScanGeometry()
        self.dataUrls = []
        with h5py.File(self.dataPath,'r') as h5In:
            for scan in self.scans:
                self.dataUrls.append(h5In.get('%s/measurement/%s'%(scan,self.xrddetector),getlink=True).path)
        
        