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

class Input:
    """
    Class to input data for processing XRD/XRF/PDF-CT and 3DXRD.
    """
    def __init__(self, dataPath):
        """
        Initialize data path.
        """
        self.data = dataPath
    
    def getScanGeometry(self):
        """
        Reads scan Geometry from h5 file.
        """
        scanKey=['1']
        scans = []
        self.geometry = []
        self.y = []
        self.rot = []
        yMotor = self.yMotor()
        rotMotor = self.rotMotor()
        with h5py.File(self.data,'r') as h5In:
            self.scans = h5In.keys()
            for scan in self.scans:
                if scan.split('.')[1] == scanKey[0]:
                    scans.append(int(scan.split('.')[0]))
        for scan in scans:
            self.geometry.append('%s.%s'%(scan,scanKey[0]))
        with h5py.File(self.data,'r') as h5In:
            for scan in self.geometry:
                self.y.append(h5In['%s/measurement/%s'%(scan,yMotor)][()])
                self.rot.append(h5In['%s/measurement/%s'%(scan,rotMotor)][()])
        print('[INFO] Read %s scans'%max(scans))
        return self.getScanGeometry
        
    def yMotor(self):
        """
        Finds y motor from provided list.
        """
        yMotors = ['dty','diffty','diffy']
        with h5py.File(self.data,'r') as h5In:
            for motor in yMotors:
                if motor in list(h5In['1.1/measurement'].keys()):
                    self.yMotor = motor
                    print('[INFO] Y motor detected: %s'%self.yMotor)
        return self.yMotor
        
    def rotMotor(self):
        """
        Finds y motor from provided list.
        """
        rotMotors = ['rot','diffrz','shrz']
        with h5py.File(self.data,'r') as h5In:
            for motor in rotMotors:
                if motor in list(h5In['1.1/measurement'].keys()):
                    self.rotMotor = motor
                    print('[INFO] Rot motor detected: %s'%self.rotMotor)
        return self.rotMotor          
                    
    def readH5(self):
        """
        Reads data from bliss HDF5.
        """
        if os.path.isfile(self.data):
            
            print(self.data)
        else:
            print('[ERROR] File not found')