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

import multiprocessing
import os

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import PyXRDCT.nmutils.utils.saveh5 as saveh5

nbprocs = int(multiprocessing.cpu_count())
try:
    nbprocs = int(os.environ['SLURM_CPUS_ON_NODE'])
except:
    print("[WARNING] Can't find SLURM_CPUS_ON_NODE")

mpl.rc('image', cmap='gray')


def shift_sino(s, s_shift):
    from scipy.ndimage import shift
    sOut = np.zeros(np.shape(s))
    for i in range(0, np.size(s, 1)):
        shift(s[:, i], -s_shift, sOut[:, i], mode='nearest')
    return sOut


def no_monitor_norm(s):
    if len(s.shape) == 3:
        avgs = np.average(s, axis=2)
        for j in range(s.shape[2]):
            for i in range(s.shape[0]):
                s[i, :, j] = s[i, :, j] / np.average(avgs[i, :])
    else:
        avgs = np.copy(s)
        for i in range(s.shape[0]):
            s[i, :] = s[i, :] / np.average(avgs[i, :])
    return (s)


class Reconstruction:
    """
    Initialise reconstruction class
    """

    def __init__(self, readH5Input, intFile=False):
        self.data = readH5Input
        if intFile:
            self.integrate = intFile

    def parallel_iradon(self, chunk):
        from skimage.transform import iradon
        recon = []
        for sino in range(chunk.shape[2]):
            recon.append(iradon(chunk[:, :, sino], sorted(self.data.rot[0]), output_size=len(self.data.y)))
        return np.array(recon)

    def parallel_histogram(self, chunk):
        import numpy as np
        sinoFinal = []
        for data in range(chunk.shape[2]):
            sino, a, y = np.histogram2d(np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),
                                        weights=np.array(chunk[:, :, data]).ravel(),
                                        bins=(self.data.rot.shape[1], self.data.rot.shape[0]))
            sinoFinal.append(sino)
        return np.array(sinoFinal)

    def reconstruct2d_s3dxrd(self, binning=1, shift=0, plot=False, save=True, no_monitor=False):
        """
        Reconstructs 2D slice of grains from segmented s3DXRD.
        """
        from skimage.transform import iradon
        tdxrdData = np.empty((len(self.data.y), len(self.data.rot[0])), dtype=np.float32)
        with h5py.File(os.path.join(self.data.savePath, 's3dxrd_segmented', self.data.dataset + '_s3dxrd_segmented.h5'),
                       'r') as h5In:
            for i, scan in enumerate(self.data.scans):
                tdxrdData[i] = h5In[scan]['nnz'][:]
        tdxrdDataSino, a, y = np.histogram2d(np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),
                                             weights=np.array(tdxrdData).ravel(), bins=(
            int(self.data.rot.shape[1] / binning), int(self.data.rot.shape[0] / binning)))
        if no_monitor:
            tdxrdDataRecon = iradon(shift_sino(no_monitor_norm(tdxrdDataSino.T), shift),
                                    sorted(self.data.rot[0][::binning]), circle=True,
                                    output_size=int(len(self.data.y) / binning))
        else:
            tdxrdDataRecon = iradon(shift_sino(tdxrdDataSino.T, shift), sorted(self.data.rot[0][::binning]),
                                    circle=True, output_size=int(len(self.data.y) / binning))
        if save:
            saveh5.saveReconstructedH5(
                os.path.join(self.data.savePath, self.data.dataset + '_s3dxrd_2dreconstruction.h5'), tdxrdDataRecon)
        if plot:
            plt.figure(figsize=(20, 10))
            plt.subplot(121)
            plt.imshow(tdxrdDataSino)
            plt.title('%s: Segmented grains sinogram' % self.data.dataset)
            plt.subplot(122)
            plt.imshow(tdxrdDataRecon)
            plt.title('%s: Segmented grains reconstruction' % self.data.dataset)
            plt.show()

    def reconstruct2d_xrdct(self, tths=[3, 4], width=0.05, binning=1, shift=0, plot=False, save=True, no_monitor=False):
        """
        Reconstructs 2D slice of XRD-CT from provided array of energies.
        """
        with h5py.File(os.path.join(self.data.savePath, 'h5_pyFAI_integrated', self.data.dataset + '_pyFAI_1.1.h5'),
                       'r') as h5In:
            tthMin = min(h5In['entry/results/polar_angle'][:])
            tthMax = max(h5In['entry/results/polar_angle'][:])
            nbptRad = len(h5In['entry/results/polar_angle'][:])
        xrdDataReconSave = []
        from skimage.transform import iradon
        for tth in tths:
            idx = (np.abs(np.linspace(tthMin, tthMax, nbptRad) - tth)).argmin()
            idxWidth = int((nbptRad / (tthMax - tthMin)) * width)
            xrdDataAvg = []
            xrdData = np.empty((len(self.data.y), len(self.data.rot[0])), dtype=np.float32)
            for i, url in enumerate(self.data.dataUrls):
                with h5py.File(os.path.join(self.data.savePath, 'h5_pyFAI_integrated',
                                            self.data.dataset + '_pyFAI_%s.h5' % (url.split('/')[1])), 'r') as h5In:
                    xrdData[i] = np.average(h5In['entry/results/data'][:, idx - idxWidth:idx + idxWidth], axis=1)
                    xrdDataAvg.append(np.average(h5In['entry/results/data'],axis=0))
            xrdDataAvg = np.average(np.array(xrdDataAvg),axis=0)
            xrdDataSino, a, y = np.histogram2d(np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),
                                               weights=np.array(xrdData).ravel(), bins=(
                int(self.data.rot.shape[1] / binning), int(self.data.rot.shape[0] / binning)))
            if no_monitor:
                xrdDataRecon = iradon(shift_sino(no_monitor_norm(xrdDataSino.T), shift),
                                      sorted(self.data.rot[0][::binning]), circle=True,
                                      output_size=int(len(self.data.y) / binning))
            else:
                xrdDataRecon = iradon(shift_sino(xrdDataSino.T, shift), sorted(self.data.rot[0]), circle=True,
                                      output_size=int(len(self.data.y) / binning))
            radius=(int(len(self.data.y) / binning)*0.9)/2
            xpr, ypr = np.mgrid[:int(len(self.data.y) / binning), :int(len(self.data.y) / binning)] - int(len(self.data.y) / binning)/2
            xrdDataReconCircle = (xpr ** 2 + ypr ** 2) > radius ** 2
            xrdDataRecon[xrdDataReconCircle] = np.average(xrdDataRecon)
            if plot:
                plt.figure(figsize=(plot, plot))
                ax1 = plt.subplot(221)
                ax1.imshow(xrdDataSino, aspect='auto')
                ax1.set_title('%s: XRD sinogram %s%s +/-%s%s' % (self.data.dataset, tth, chr(176), width, chr(176)))
                ax2 = plt.subplot(222)
                ax2.imshow(xrdDataRecon,vmin=1.2*np.min(xrdDataRecon),vmax=0.8*np.max(xrdDataRecon))
                ax2.set_title('%s: XRD reconstruction %s%s +/-%s%s' % (self.data.dataset, tth, chr(176), width, chr(176)))
                ax3 = plt.subplot(212)
                ax3.plot(np.linspace(tthMin, tthMax, nbptRad),np.log(xrdDataAvg),'k')
                ax3.set_xlabel('tth (%s)'%chr(176))
                ax3.set_ylabel('log(I) (A.U.)')
                ax3.axvline(x = tth-width, color = 'r')
                ax3.axvline(x = tth+width, color = 'r')
                ax3.set_title('%s: XRD pattern (average)'% (self.data.dataset))
                plt.show()
            with open(os.path.join(self.data.savePath, self.data.dataset + '_xrd_sum.xy'), "w" ) as f:
                f.write("#  tth(deg)  intensity\n")
                for twoth,intens in zip(np.linspace(tthMin, tthMax, nbptRad),xrdDataAvg):
                    f.write("%g %g\n"%(twoth,intens))
            xrdDataReconSave.append(xrdDataRecon)
        xrdDataReconSave = np.array(xrdDataReconSave)
        if save:
            saveh5.saveReconstructedH5(os.path.join(self.data.savePath, self.data.dataset + '_xrd_2dreconstruction.h5'),
                                       xrdDataReconSave, tths, xAxis='tth')

    def reconstruct3d_xrdct(self, algorithm='fbp', binning=1, shift=0, save=True, no_monitor=False,plot=False):
        """
        Reconstructs 3D dataset of XRD-CT from provided array of energies.
        """
        if not os.path.exists(os.path.join(self.data.savePath, self.data.dataset + '_xrd_3dreconstruction.h5')):
            with h5py.File(os.path.join(self.data.savePath, 'h5_pyFAI_integrated', self.data.dataset + '_pyFAI_1.1.h5'),
                           'r') as h5In:
                tth = h5In['entry/results/polar_angle'][:]
            xrdData = np.empty((len(self.data.y), len(self.data.rot[0]), len(tth)), dtype=np.float32)
            for i, url in enumerate(self.data.dataUrls):
                with h5py.File(os.path.join(self.data.savePath, 'h5_pyFAI_integrated',
                                            self.data.dataset + '_pyFAI_%s.h5' % (url.split('/')[1])), 'r') as h5In:
                    xrdData[i, :, :] = h5In['entry/results/data'][:]
            xrdDataSino = np.zeros(
                (int(self.data.rot.shape[0] / binning), int(self.data.rot.shape[1] / binning), xrdData.shape[2]))
            for tthVal in range(xrdData.shape[2]):
                xrdDataSinoBuffer, a, y = np.histogram2d(np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),
                                                         weights=np.array(xrdData[:, :, tthVal]).ravel(), bins=(
                    int(self.data.rot.shape[1] / binning), int(self.data.rot.shape[0] / binning)))
                xrdDataSino[:, :, tthVal] = shift_sino(xrdDataSinoBuffer.T, shift)
            chunks = [xrdDataSino[:, :, proc * 100:(proc * 100) + 100] for proc in range(int(xrdData.shape[2] / 100))]
            xrdDataReconSave = []
            if no_monitor:
                xrdDataSino = no_monitor_norm(xrdDataSino)
            with multiprocessing.Pool(int(multiprocessing.cpu_count() / 2)) as pool:
                for result in pool.map(self.parallel_iradon, chunks):
                    xrdDataReconSave.extend(result)
            xrdDataReconSave = np.array(xrdDataReconSave)
            xrdDataSino = xrdDataSino.T
            if save:
                saveh5.saveReconstructedH5(os.path.join(self.data.savePath, self.data.dataset + '_xrd_3dreconstruction.h5'),
                                           xrdDataReconSave, tth, xAxis='tth')
                saveh5.saveReconstructedH5(os.path.join(self.data.savePath, self.data.dataset + '_xrd_3dsinogram.h5'),
                                           xrdDataSino, tth, xAxis='tth')
        else:
            print('[INFO] Found already reconstructed datasets!')
            with h5py.File(os.path.join(self.data.savePath, self.data.dataset + '_xrd_3dreconstruction.h5'),'r') as h5In:
                xrdDataReconSave = h5In['entry_0000/data'][:]
            with h5py.File(os.path.join(self.data.savePath, self.data.dataset + '_xrd_3dsinogram.h5'),'r') as h5In:
                xrdDataSino = h5In['entry_0000/data'][:]
                tth = h5In['entry_0000/tth'][:]
        if plot:
            xrdDataAvg = np.average(np.average(xrdDataSino,axis=1),axis=1)
            from matplotlib.widgets import Slider, Button
            def update_tth(val):
                idx = (np.abs(np.linspace(tth[0], tth[-1], len(tth)) - tth_slider.val)).argmin()
                ax1.imshow(xrdDataSino[idx,:,:], aspect='auto')
                ax1.set_title('%s: XRD sinogram %.2f%s' % (self.data.dataset, tth[idx], chr(176)))
                ax2.imshow(xrdDataReconSave[idx,:,:])#,vmin=1.2*np.min(xrdDataReconSave),vmax=0.8*np.max(xrdDataReconSave))
                ax2.set_title('%s: XRD reconstruction %.2f%s' % (self.data.dataset, tth[idx], chr(176)))
                ax3.set_ydata(tth_slider.val)
                fig.canvas.draw_idle()
            fig = plt.figure(figsize=(plot, plot*0.66))
            fig.subplots_adjust(bottom=0.22)
            ax1 = plt.subplot(221)
            ax1.imshow(xrdDataSino[100,:,:], aspect='auto')
            ax1.set_title('%s: XRD sinogram %.2f%s' % (self.data.dataset, tth[0], chr(176)))
            ax2 = plt.subplot(222)
            ax2.imshow(xrdDataReconSave[100,:,:])#,vmin=1.2*np.min(xrdDataReconSave),vmax=0.8*np.max(xrdDataReconSave))
            ax2.set_title('%s: XRD reconstruction %.2f%s' % (self.data.dataset, tth[0], chr(176)))
            ax3 = plt.subplot(212)
            ax3.plot(np.linspace(tth[0], tth[-1], len(tth)),np.log(xrdDataAvg),'k')
            ax3.set_xlabel('2\u03B8 (%s)'%chr(176))
            ax3.set_ylabel('log(I) (A.U.)')
            #ax4 = fig.add_axes([0.15, 0.1, 0.8, 0.03])
            #handle_style = dict(style='x')
            tth_slider = Slider(ax=ax3, label='2\u03B8 (%s)'%chr(176), valmin=float(tth[0]), valmax = float(tth[-1]),valinit=float(tth[100]),orientation="horizontal",color='grey')
            tth_slider.on_changed(update_tth)
            plt.show()
            
    def reconstruct2d_xrfct(self, energies=[2.013, 28.612, 49.127, 61.140, 0.5249, 0.0543, 4.952, 28.612, 49.127],
                            binning=1, width=0.05, shift=0, plot=False, save=True, no_monitor=False):
        """
        Reconstructs 2D slice of XRF-CT from provided array of energies.
        """
        from skimage.transform import iradon
        ENERGY_MAX = 81.92
        xrfDataReconSave = []
        for energy in energies:
            idx = (np.abs(np.linspace(0, ENERGY_MAX, self.data.channels) - energy)).argmin()
            idxWidth = int((self.data.channels / ENERGY_MAX) * width)
            xrfData = np.empty((len(self.data.y), len(self.data.rot[0])), dtype=np.float32)
            with h5py.File(self.data.dataPath, 'r') as h5In:
                for i, scan in enumerate(self.data.scans):
                    xrfData[i] = np.average(
                        h5In[scan]['measurement'][self.data.xrfdetector][:, idx - idxWidth:idx + idxWidth], axis=1) / \
                                 h5In[scan]['measurement'][self.data.beamMonitor][:]
            xrfDataSino, a, y = np.histogram2d(np.array(self.data.rot).ravel(), np.array(self.data.y).ravel(),
                                               weights=np.array(xrfData).ravel(), bins=(
                int(self.data.rot.shape[1] / binning), int(self.data.rot.shape[0] / binning)))
            if no_monitor:
                xrfDataRecon = iradon(shift_sino(no_monitor_norm(xrfDataSino.T), shift),
                                      sorted(self.data.rot[0][::binning]), circle=True,
                                      output_size=int(len(self.data.y) / binning))
            else:
                xrfDataRecon = iradon(shift_sino(xrfDataSino.T, shift), sorted(self.data.rot[0]), circle=True,
                                      output_size=int(len(self.data.y) / binning))
            if plot:
                plt.figure(figsize=(20, 10))
                plt.subplot(121)
                plt.imshow(xrfDataSino)
                plt.title('%s: XRF sinogram %skeV +/-%skeV' % (self.data.dataset, energy, width))
                plt.subplot(122)
                plt.imshow(xrfDataRecon)
                plt.title('%s: XRF reconstruction %skeV +/-%skeV' % (self.data.dataset, energy, width))
                plt.show()
            xrfDataReconSave.append(xrfDataRecon)
        xrfDataReconSave = np.array(xrfDataReconSave)
        if save:
            saveh5.saveReconstructedH5(os.path.join(self.data.savePath, self.data.dataset + '_xrf_2dreconstruction.h5'),
                                       xrfDataReconSave, energies, xAxis='Energy')

    def reconstruct3d_xrfct(self, algorithm='fbp', binning=1, shift=0, save=True, no_monitor=False):
        """
        Reconstructs 3D dataset of XRF-CT from provided array of energies.
        """
        import multiprocessing
        ENERGY_MAX = 81.92
        energies = np.linspace(0, ENERGY_MAX, self.data.channels)
        xrfData = np.empty((len(self.data.y), len(self.data.rot[0]), len(energies)), dtype=np.float32)
        with h5py.File(self.data.dataPath, 'r') as h5In:
            for i, scan in enumerate(self.data.scans):
                xrfData[i] = h5In[scan]['measurement'][self.data.xrfdetector][:] / h5In[scan]['measurement'][
                                                                                       self.data.beamMonitor][:][:,
                                                                                   None]
        dataChunks = [xrfData[:, :, proc * 100:(proc * 100) + 100] for proc in range(int(xrfData.shape[2] / 100))]
        xrfDataSino = []
        with multiprocessing.Pool(nbprocs) as pool:
            for result in pool.map(self.parallel_histogram, dataChunks):
                xrfDataSino.extend(result)
        xrfDataSino = np.array(xrfDataSino)
        for energyVal in range(xrfData.shape[2]):
            xrfDataSino[:, :, energyVal] = shift_sino(xrfDataSino.T, shift)
        sinoChunks = [xrfDataSino[:, :, proc * 100:(proc * 100) + 100] for proc in range(int(xrfData.shape[2] / 100))]
        xrfDataReconSave = []
        if no_monitor:
            xrfDataSino = no_monitor_norm(xrfDataSino)
        with multiprocessing.Pool(nbprocs) as pool:
            for result in pool.map(self.parallel_iradon, sinoChunks):
                xrfDataReconSave.extend(result)
        xrfDataReconSave = np.array(xrfDataReconSave)
        if save:
            saveh5.saveReconstructedH5(os.path.join(self.data.savePath, self.data.dataset + '_xrf_3dreconstruction.h5'),
                                       xrfDataReconSave, energies, xAxis='energy')
            saveh5.saveReconstructedH5(os.path.join(self.data.savePath, self.data.dataset + '_xrf_3dsinogram.h5'),
                                       xrfDataSino, energies, xAxis='energy')
