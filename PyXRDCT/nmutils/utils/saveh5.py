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

import os

import h5py
from pyFAI.io.nexus import save_NXmonpd


def makeSaveDirs(savePath):
    if not os.path.exists(savePath):
        os.makedirs(savePath)
        print('[INFO] Folder %s created' % savePath)


def saveIntegrateH5(savePath, result, title):
    """
    Saves result in a h5 NeXus file with relpath being filled with the rest of the save path (including extension)
    """
    makeSaveDirs(os.path.dirname(savePath))
    save_NXmonpd(savePath, result, title=title, entry='entry', instrument='ID11 beamline', source_name='ESRF',
                 source_type='Synchrotron', source_probe='x-ray',
                 sample=os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(savePath)))), extra=None)


def saveReconstructedH5(savePath, result, metadata=None, xAxis='X'):
    """
    Saves reconstruction data with metadata as h5
    """
    if metadata is None:
        metadata = []
    makeSaveDirs(os.path.dirname(savePath))
    with h5py.File(savePath, 'w') as h5Out:
        dsetResult = h5Out.create_dataset('entry_0000/data', result.shape, dtype='f')
        dsetResult[...] = result
        dsetMetadata = h5Out.create_dataset('entry_0000/%s' % xAxis, [len(metadata)], dtype='f')
        dsetMetadata[...] = metadata
    print('[INFO] %s saved!' % savePath)
