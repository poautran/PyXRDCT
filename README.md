# PyXRDCT: X-Ray Diffraction Computed Tomography reconstruction functions

PyXRDCT provides X-Ray Diffraction Computed Tomography reconstruction functions.
The library gives functions for correcting diffraction tomography sinograms based on simple mathematics such as center of mass or outliers values.

## Installation

### Cloning

The source code is available here by running the command:

	git clone https://github.com/poautran/PyXRDCT.git

Then run the following commands to install it:

	cd PyXRDCT
	python setup.py sdist bdist_wheel

The package should be working, to test it run:

	ipython
	import PyXRDCT

## Dependencies and OS

The library is currently tested on Python 3.6.8 over Ubuntu 18.04. Further version will be tested in the future.

### Requirements :

	* numpy 		- 	http://www.numpy.org
	* scipy 		- 	http://www.scipy.org
	* matplotlib 		- 	http://matplotlib.sourceforge.net/
	* h5py	    		-  	http://www.h5py.org/
	* fabio			-	https://github.com/silx-kit/fabio
	* pyFAI			-	https://pyfai.readthedocs.io/
	* Others, not detailed yet

## Documentation

The user manual is under construction. For an early start, main code to be run is in PyXRDCT/core. Functions and other utilities are in PyXRDCT/nmutils/utils. Examples of raw data are located in the resources folder in PyXRDCT/resources.

To start you can add a shortcut on to your .bashrc file such as:

	alias pyxrdct='/path/to/folder/PyXRDCT/PyXRDCT/core.py'

Then run the following command to see all the options:

	pyxrdct -h



