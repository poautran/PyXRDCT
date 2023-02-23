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

import json
import fabio

import multiprocessing

def integrator(urls, config, data):
    import os, fabio, pyFAI, pyFAI.azimuthalIntegrator as AI, numpy as np, h5py, hdf5plugin, PyXRDCT.PyXRDCT.nmutils.utils.saveh5 as saveh5
    for url in urls:
        saveIntH5Path = os.path.join(data.savePath,'h5_pyFAI_integrated',data.dataset+'_pyFAI_%s.h5'%(url.split('/')[1]))
        mask = fabio.open(config['mask_file']).data
        method = pyFAI.method_registry.IntegrationMethod.select_method(dim=1, split="full", algo="csr", impl="opencl")[0]
        resultFinal = []
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
        with h5py.File(data.dataPath,'r') as h5In:
            result = np.empty((h5In[url].shape[0], config['nbpt_rad']), dtype=np.float32)
            monitor = h5In[url.split('/')[1]]['measurement'][data.beamMonitor][:]*1e-6
            for image in range(h5In[url].shape[0]):
                resultBuffer = ai.integrate1d_ng(h5In[url][image,:,:],
                                      config['nbpt_rad'],
                                      mask=mask,
                                      method=method,
                                      unit=config['unit'],
                                      polarization_factor=config['polarization_factor'],
                                      correctSolidAngle=config['do_solid_angle']
                                      )
                result[image,:] = resultBuffer.intensity/monitor[image]
        radial = resultBuffer.radial
        saveh5.saveIntegrateH5(saveIntH5Path,resultBuffer,'XRDCT: pyFAI integration scan %s'%(url.split('/')[1]))
        print('[INFO] %s DONE!'%saveIntH5Path)
        with h5py.File(saveIntH5Path, 'r+') as h5In:
            del h5In['entry/results/data']
            h5In.create_dataset('entry/results/data',data=result)
        
#def integrator_bitshuffle(urls, config, data):
#    import os, fabio, pyFAI, pyFAI.azimuthalIntegrator as AI, numpy as np, h5py, hdf5plugin, saveh5
#    for url in urls:
#        saveIntH5Path = os.path.join(data.savePath,'h5_pyFAI_integrated',data.dataset+'_pyFAI_%s.h5'%(url.split('/')[1]))
#        mask = fabio.open(config['mask_file']).data
#        method = pyFAI.method_registry.IntegrationMethod.select_method(dim=1, split="full", algo="csr", impl="opencl")[0]
#        resultFinal = []
#        ai = AI.AzimuthalIntegrator(dist=config['dist'], 
#                                    poni1=config['poni1'], 
#                                    poni2=config['poni2'], 
#                                    rot1=config['rot1'], 
#                                    rot2=config['rot2'], 
#                                    rot3=config['rot3'], 
#                                    detector=pyFAI.detectors.NexusDetector(config['detector_config']['filename']), 
#                                    wavelength=config['wavelength']) 
#        with h5py.File(os.path.join(os.path.dirname(data.dataPath),'scan%04d/%s_0000.h5'%(float(url.split('/')[1]),data.xrddetector)),'r') as h5In:
#            ds = h5In["entry_0000/measurement/data"]
#            print(ds.id)
#            filter_mask, chunk = ds.id.read_direct_chunk(ds.id.get_chunk_info(0).chunk_offset)
#            gpu_decompressor = BitshuffleLz4(len(chunk), frame.size, dtype=frame.dtype, ctx=engine.ctx)
#            gpu_decompressor.block_size = 128 
#            result = np.empty((h5In[url].shape[0], config['nbpt_rad']), dtype=np.float32)
#            for image in range(h5In[url].shape[0]):
#                resultBuffer = engine.integrate_ng(gpu_decompressor(chunk),
#                                      config['nbpt_rad'],
#                                      mask=mask,
#                                      method=method,
#                                      unit=config['unit'],
#                                      polarization_factor=config['polarization_factor'],
#                                      correctSolidAngle=config['do_solid_angle']
#                                      )
#                result[image,:] = resultBuffer.intensity
#        radial = resultBuffer.radial
#        saveh5.saveIntegrateH5(saveIntH5Path,resultBuffer,'XRDCT: pyFAI integration scan %s'%(url.split('/')[1]))
#        print('[INFO] %s DONE!'%saveIntH5Path)
#        with h5py.File(saveIntH5Path, 'r+') as h5In:
#            del h5In['entry/results/data']
#            h5In.create_dataset('entry/results/data',data=result)
        
class Integrate:
    """
    Initialise integrate class
    """

    def __init__(self, readH5Input, jsonFile):
        self.data = readH5Input
        with open(jsonFile) as jsonIn:
            self.config = json.load(jsonIn)    

    def wrap(self, chunk ):
        integrator( chunk, self.config, self.data)
    
    def integrate1d(self):
        chunks = [self.data.dataUrls[proc::int(multiprocessing.cpu_count()/2)] for proc in range(int(multiprocessing.cpu_count()/2))]
        start_time = time.time()
        with multiprocessing.Pool(int(multiprocessing.cpu_count()/2)) as pool:
            for _ in pool.imap_unordered( self.wrap, chunks ):
                pass
        print('[INFO] Took: %4dsec, %4dFPS'%(time.time()-start_time,(len(self.data.dataUrls)*len(self.data.rot[0,:]))/(time.time()-start_time)))

    
