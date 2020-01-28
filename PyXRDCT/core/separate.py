#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, glob
import matplotlib.pyplot as plt
import h5py
import numpy as np
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import argparse
import re
import fabio

### Functions ###

def readJson(jsonFile):
    import json
    with open(jsonFile) as f:
        config = json.load(f)
    return config


### Input ###

def main():
	parser = argparse.ArgumentParser(description='This script separates the powder from monocristalline part on 2D images. python separate.py FILE -j azimint.json')
	parser.add_argument('INPUT',metavar="INPUT", help="List of files. Ex : *.edf", nargs='+')
	parser.add_argument('-j','--json',dest="JSON", default=".azimint.json", help="Configuration file containing the processing to be done")
	parser.add_argument('-o','--output',type=str,dest='OUTPUT',help='Output path')
	parser.add_argument('-d','--display',help='Display result images and integrated patterns',dest='DISPLAY',action='store_true')
	parser.set_defaults(func=run)
	args = parser.parse_args()
	args.func(args)


def run(args):
	FILE = args.INPUT
	if np.shape(FILE)==(1,):
		FILE = np.sort(glob.glob(str(FILE[0])))

	SAVE_PATH = os.getcwd()
	if args.OUTPUT:
		SAVE_PATH = args.OUTPUT
	else:
		print('!!! Warning files will be saved in the current folder because no output was defined.')

	### Read json file ###
	jsonParam = readJson(args.JSON)
	### Integration of FILE ###
	azimutalIntegrator = AzimuthalIntegrator(dist=jsonParam['dist'], poni1=jsonParam['poni1'], poni2=jsonParam['poni2'], rot1=jsonParam['rot1'], rot2=jsonParam['rot2'], rot3=jsonParam['rot3'], pixel1=jsonParam['detector_config']['pixel1'], pixel2=jsonParam['detector_config']['pixel2'], detector=jsonParam['detector'], wavelength=jsonParam['wavelength'])
	#dark = np.array(fabio.open(jsonParam['dark_current']).data)
	#flat = np.array(fabio.open(jsonParam['flat_field']).data)
	mask = np.array(fabio.open(jsonParam['mask_file']).data)
	for i in range(len(FILE)):
		dataFile = np.array(fabio.open(FILE[i]).data)
		dataBragg, dataPowder = azimutalIntegrator.separate(dataFile, npt_rad=1024, npt_azim=512, unit=jsonParam['unit'], method='splitpixel', percentile=50, mask=None, restore_mask=True)
		dataXP, dataYP = azimutalIntegrator.integrate1d(dataPowder, int(jsonParam['nbpt_rad']), filename=None, correctSolidAngle=jsonParam['do_solid_angle'], variance=None, error_model=None, radial_range=(float(jsonParam['radial_range_min']),float(jsonParam['radial_range_max'])), azimuth_range=None, mask=mask, dummy=jsonParam['do_dummy'], delta_dummy=jsonParam['delta_dummy'], polarization_factor=None, method='csr', unit=jsonParam['unit'], safe=True, profile=False, all=False, metadata=None)
		dataXB, dataYB = azimutalIntegrator.integrate1d(dataBragg, int(jsonParam['nbpt_rad']), filename=None, correctSolidAngle=jsonParam['do_solid_angle'], variance=None, error_model=None, radial_range=(float(jsonParam['radial_range_min']),float(jsonParam['radial_range_max'])), azimuth_range=None, mask=mask, dummy=jsonParam['do_dummy'], delta_dummy=jsonParam['delta_dummy'], polarization_factor=None, method='csr', unit=jsonParam['unit'], safe=True, profile=False, all=False, metadata=None)
		dataXO, dataYO = azimutalIntegrator.integrate1d(dataFile, int(jsonParam['nbpt_rad']), filename=None, correctSolidAngle=jsonParam['do_solid_angle'], variance=None, error_model=None, radial_range=(float(jsonParam['radial_range_min']),float(jsonParam['radial_range_max'])), azimuth_range=None, mask=mask, dummy=jsonParam['do_dummy'], delta_dummy=jsonParam['delta_dummy'], polarization_factor=None, method='csr', unit=jsonParam['unit'], safe=True, profile=False, all=False, metadata=None)
		
		# Display #
		if args.DISPLAY:
			ax1 = plt.subplot(2, 3, 1)
			ax1.imshow(dataFile)
			ax1.set_title('Original data')	
			ax2 = plt.subplot(2, 3, 2)	
			ax2.imshow(dataBragg)
			ax2.set_title('Bragg contribution')	
			ax3 = plt.subplot(2, 3, 3)
			ax3.imshow(dataPowder)
			ax3.set_title('Powder contribution')	
			ax4 = plt.subplot(2, 1, 2)	
			ax4.plot(dataXP, dataYP,label='Powder contribution')
			ax4.plot(dataXB, dataYB,label='Bragg contribution')
			ax4.plot(dataXO, dataYO,label='Original')
			ax4.set_title('Integrated data')
			ax4.set_xlabel(r'$Q (\AA^{-1})$', fontsize=11)
			ax4.set_ylabel(r'$I(A.U.)$', fontsize=11)
			plt.legend()
			plt.show()

		# Save Edf File #		

		#f.open(FILE[:-4]+'powder_contrib.edf')




		# Save integrated data #
		currentFILE = FILE[i]
		f_out_dataPowder = open(currentFILE[:-4] + '_powder.dat','w')
		for j in range(0,np.size(dataXP,0)):
			f_out_dataPowder.writelines('%07f'%dataXP[j] + '    ' + '%04f'%dataYP[j]  + '\n')
		f_out_dataPowder.close()
		f_out_dataBragg = open(currentFILE[:-4] + '_bragg.dat','w')
		for j in range(0,np.size(dataXB,0)):
			f_out_dataBragg.writelines('%07f'%dataXB[j] + '    ' + '%04f'%dataYB[j]  + '\n')
		f_out_dataBragg.close()
		f_out_dataOriginal = open(currentFILE[:-4] + '_original.dat','w')
		for j in range(0,np.size(dataXO,0)):
			f_out_dataOriginal.writelines('%07f'%dataXO[j] + '    ' + '%04f'%dataYO[j]  + '\n')
		f_out_dataOriginal.close()



if __name__=="__main__":
	main()


