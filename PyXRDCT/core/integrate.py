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
# FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import json
import multiprocessing
import time
import os

nbprocs = int(multiprocessing.cpu_count())
try:
    nbprocs = int(os.environ['SLURM_CPUS_ON_NODE'])
except:
    print("[WARNING] Can't find SLURM_CPUS_ON_NODE")

def integrator(urls, jsonPath, data):
    global detector
    import os, fabio, json, multiprocessing, pyFAI, pyFAI.azimuthalIntegrator as AI, numpy as np, time, hdf5plugin, h5py, \
        PyXRDCT.nmutils.utils.saveh5 as saveh5
    os.environ["OMP_NUM_THREADS"] = "1"
    with open(jsonPath) as jsonIn:
        config = json.load(jsonIn)
    startTime = time.time()
    mask = fabio.open(config['mask_file']).data
    if config['do_dark']:
        dark = fabio.open(config['dark_current'][0]).data
    else:
        dark = None
    if config['do_flat']:
        flat = fabio.open(config['flat_field'][0]).data
    else:
        flat = None
    if config['do_radial_range']:
        radial_range = (config['radial_range_min'], config['radial_range_max'])
    else:
        radial_range = None
    if config['do_azimuthal_range']:
        azimuth_range = (config['azimuth_range_min'], config['azimuth_range_max'])
    else:
        azimuth_range = None
    method = pyFAI.method_registry.IntegrationMethod.select_method(dim=1, split="full", algo="csc", impl="cython")[0]
    if 'filename' in config['detector_config'].keys():
        detector = pyFAI.detectors.NexusDetector(config['detector_config']['filename'])
    elif 'splineFile' in config['detector_config'].keys():
        detector = pyFAI.detectors.Detector(splineFile=config['detector_config']['splineFile'])
    else:
        print('[WARNING] Detector config not found!')
    ai = AI.AzimuthalIntegrator(dist=config['dist'],
                                poni1=config['poni1'],
                                poni2=config['poni2'],
                                rot1=config['rot1'],
                                rot2=config['rot2'],
                                rot3=config['rot3'],
                                detector=detector,
                                wavelength=config['wavelength'])
    for url in urls:
        os.sched_setaffinity(0, [int(float(int(url.split('/')[1].split('.')[0])-1) % multiprocessing.cpu_count())])
        saveIntH5Path = os.path.join(data.savePath, 'h5_pyFAI_integrated', data.dataset + '_pyFAI_%s.h5' % (url.split('/')[1]))
        checkDoneIntegration = os.path.exists(saveIntH5Path)
        if checkDoneIntegration:
            print('%s Already processed!'%url)
            continue
        with h5py.File(data.dataPath, 'r') as h5In:
            frameStackShape = [h5In[url].shape[0],h5In[url].shape[1],h5In[url].shape[2]]
            result = np.empty((h5In[url].shape[0], config['nbpt_rad']), dtype=np.float32)
            monitor = h5In[url.split('/')[1]]['measurement'][data.beamMonitor][:] * 1e-6
        resultBuffer = []
        readBuffer = np.zeros((frameStackShape[1],frameStackShape[2]),dtype='uint32')
        with h5py.File(os.path.join(os.path.dirname(data.dataPath),'scan%04d/%s_0000.h5'%(int(url.split('/')[1].split('.')[0]),data.xrddetector)), 'r') as h5In:
            for image in range(frameStackShape[0]):
                h5In['entry_0000/measurement/data'].read_direct(readBuffer, np.s_[image,:,:], np.s_[:,:])
                resultBuffer.append(ai.integrate1d_ng(readBuffer,
                                                 config['nbpt_rad'],
                                                 mask=mask,
                                                 method=method,
                                                 dark=dark,
                                                 flat=flat,
                                                 radial_range=radial_range,
                                                 azimuth_range=azimuth_range,
                                                 polarization_factor=float(config['polarization_factor']),
                                                 unit=config['unit']
                                                 ).intensity
                                   )
            result = resultBuffer / monitor[:,None]
        resultSave = ai.integrate1d_ng(readBuffer,config['nbpt_rad'],mask=mask,method=method,dark=dark,flat=flat,radial_range=radial_range,azimuth_range=azimuth_range,polarization_factor=float(config['polarization_factor']),unit=config['unit'])
        saveh5.saveIntegrateH5(saveIntH5Path, resultSave, 'XRDCT: pyFAI integration scan %s' % (url.split('/')[1]))
        print('[INFO] %s DONE! Took %s seconds!' %(saveIntH5Path,time.time()-startTime))
        with h5py.File(saveIntH5Path, 'r+') as h5In:
            del h5In['entry/results/data']
            h5In.create_dataset('entry/results/data', data=result)

class Integrate:
    """
    Initialise integrate class
    """

    def __init__(self, readH5Input, jsonFile):
        self.data = readH5Input
        self.jsonPath = jsonFile
        with open(jsonFile) as jsonIn:
            self.config = json.load(jsonIn)

    def wrap(self, chunk):
        integrator(chunk, self.jsonPath, self.data)

    def integrate1d(self):
        chunks = [self.data.dataUrls[proc::nbprocs] for proc in
                  range(nbprocs)]
        start_time = time.time()
        with multiprocessing.Pool(nbprocs) as pool:
            for _ in pool.imap_unordered(self.wrap, chunks):
                pass
        print('[INFO] Took: %4dsec, %4dFPS' % (
            time.time() - start_time,
            (len(self.data.dataUrls) * len(self.data.rot[0, :])) / (time.time() - start_time)))

            

        
        
        
        
        
        
        
        
        
        
        
        
        
        