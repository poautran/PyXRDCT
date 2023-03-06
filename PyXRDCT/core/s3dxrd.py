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

#This portion of the package is based on ImageD11 from Jonathan Wright: https://github.com/FABLE-3DXRD/ImageD11

import numpy as np
import h5py, hdf5plugin
import os, sys, time

from ImageD11 import sparseframe, cImageD11
import numba
import fabio
import concurrent.futures
import multiprocessing

msk = 1 - fabio.open('/data/id11/nanoscope/Eiger/mask_20210428.edf').data  
howmany =  10000
pixels_in_spot = 5
thresholds = (4,8,16,32,64,128,256) 
CUT = 25


class bgsub( object ):
    def __init__(self, gain=0.5, sigmap=0.2, sigmat=15 ):
        self.sigmap = sigmap
        self.sigmat = sigmat
        self.gain = gain
        self.bg = None
    def __call__(self, input):
        if self.bg is None:
            self.msk = np.empty( input.shape, np.uint8 )
            self.bg = np.empty( input.shape, np.float32 )
        cImageD11.bgcalc( input,
                          self.bg,
                          self.msk,
                          self.gain,
                          self.sigmap,
                          self.sigmat )
        return input - self.bg

@numba.njit
def select( img, msk, row, col, val ):
    cut = CUT
    # Choose the pixels that are > cut and put into sparse arrays
    k = 0
    for s in range(img.shape[0]):
        for f in range(img.shape[1]):
            if img[s,f]>cut:
                if msk[s,f]: # skip masked
                    continue
                row[k] = s
                col[k] = f
                val[k] = img[s,f]
                k += 1
    return k

@numba.njit
def top_pixels( nnz, row, col, val, howmany ):
    """
    selects the strongest pixels from a sparse collection
    - thresholds should be a sorted array of potential cutoff values to try
    that are higher than the original cutoff used to select data
    - howmany is the maximum number of pixels to return
    """
    # quick return if there are already few enough pixels
    global thresholds
    if nnz <= howmany:
        return nnz
    # histogram of how many pixels are above each threshold
    h = np.zeros(len(thresholds), dtype=np.uint32)
    for k in range(nnz):
        for i,t in enumerate(thresholds):
            if val[k] > t:
                h[i] += 1
            else:
                break
    # choose the one to use. This is the first that is lower than howmany
    tcut = thresholds[-1]
    for n,t in zip( h, thresholds ):
        if n < howmany:
            tcut = t
            break
    # now we filter the pixels
    n = 0
    for k in range(nnz):
        if val[k] > tcut:
            row[n] = row[k]
            col[n] = col[k]
            val[n] = val[k]
            n += 1
            if n >= howmany:
                break
    return n

def choose_parallel( args ):
    """ reads a frame and sends back a sparse frame """
    h5name, address, frame_num = args
    with h5py.File( h5name, "r" ) as h:
        frm = h[address][frame_num]
    row = np.empty(msk.size, np.uint16)
    col = np.empty(msk.size, np.uint16)
    val = np.empty(msk.size, frm.dtype)
    nnz = select( frm, msk, row, col, val)
    global howmany
    if nnz == 0:
        sf = None
    else:
        if nnz > howmany:
            nnz = top_pixels( nnz, row, col, val, howmany )
        # Now get rid of the single pixel 'peaks'
        s = sparseframe.sparse_frame( row[:nnz].copy(), col[:nnz].copy(),
                                      frm.shape )
        s.set_pixels("intensity", val[:nnz].copy())
        # label them according to the connected objects
        sparseframe.sparse_connected_pixels( s, threshold = 5,
                                                 data_name='intensity',
                                                 label_name='cp')
        # only keep spots with more than 3 pixels ...
        mom = sparseframe.sparse_moments( s, intensity_name='intensity',
                                                 labels_name='cp' )
        npx = mom[:,cImageD11.s2D_1]
        pxcounts = npx[s.pixels['cp']-1]
        pxmsk = pxcounts >= pixels_in_spot
        if pxmsk.sum() == 0:
            sf = None
        else:
            sf = s.mask( pxmsk )
    return frame_num, sf

def segment_scans(h5FileIn):
    """ Does segmentation on a series of scans in hdf files:
    """
    opts = { 'chunks' : (10000,), 'maxshape' : (None,),'compression':'lzf', 'shuffle': True }
    ndone = 0
    outname = os.path.join(h5FileIn.savePath,'s3dxrd_segmented',h5FileIn.dataset+'_s3dxrd_segmented.h5')
    if not os.path.exists(os.path.dirname(outname)):
        os.makedirs(os.path.dirname(outname))
    with h5py.File( outname, "w") as hout:
        for scan in h5FileIn.scans:
            if scan.endswith(".2"): # for fscans
                continue
            with h5py.File(h5FileIn.dataPath, "r") as hin:
                hout.attrs['h5input'] = h5FileIn.dataPath
                hdr = {}
                g = hout.create_group( scan )
                gm = g.create_group('measurement')
                gm.create_dataset(h5FileIn.rotMotor, data = hin[scan]['measurement'][h5FileIn.rotMotor][:] )
                gm.create_dataset(h5FileIn.beamMonitor, data = hin[scan]['measurement'][h5FileIn.beamMonitor][:] )
                gip = g.create_group('instrument/positioners')
                gip.create_dataset(h5FileIn.yMotor, data = hin[scan]['instrument/positioners'][h5FileIn.yMotor][()] )
                frms = hin[scan]['measurement'][h5FileIn.xrddetector]
                row = g.create_dataset( 'row', (1,), dtype=np.uint16, **opts )
                col = g.create_dataset( 'col', (1,), dtype=np.uint16, **opts )
                # can go over 65535 frames in a scan 
                num = g.create_dataset( 'frame', (1,), dtype=np.uint32,  **opts )
                sig = g.create_dataset( 'intensity', (1,), dtype=frms.dtype, **opts )
                nnz = g.create_dataset( 'nnz', (frms.shape[0],), dtype=np.uint32)
                g.attrs['itype']  = np.dtype(np.uint16).name
                g.attrs['nframes']= frms.shape[0]
                g.attrs['shape0'] = frms.shape[1]
                g.attrs['shape1'] = frms.shape[2]
                npx = 0               
                address = scan+"/measurement/"+h5FileIn.xrddetector
                nimg = frms.shape[0]
                args = [ (h5FileIn.dataPath, address, i) for i in range(nimg) ]
            chunksize = max(1,len(args)//multiprocessing.cpu_count()//8)
            for i, spf in concurrent.futures.ProcessPoolExecutor(max_workers = multiprocessing.cpu_count()).map(choose_parallel, args, chunksize=chunksize, timeout=60):
                if spf is None:
                    nnz[i] = 0
                    continue
                if spf.nnz + npx > len(row):
                    row.resize( spf.nnz + npx, axis=0 )
                    col.resize( spf.nnz + npx, axis=0 )
                    sig.resize( spf.nnz + npx, axis=0 )
                    num.resize( spf.nnz + npx, axis=0 )
                row[npx:] = spf.row[:]
                col[npx:] = spf.col[:]
                sig[npx:] = spf.pixels['intensity']
                num[npx:] = i
                nnz[i] = spf.nnz
                npx += spf.nnz
            ndone += nimg
            try:
                print("[INFO] Scan %s DONE! Found %d spots"%(scan, spf.nnz))
            except:
                print("[INFO] Scan %s DONE! Found 0 spot"%(scan))
            sys.stdout.flush()
    return ndone
