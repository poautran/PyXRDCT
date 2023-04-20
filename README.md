# PyXRDCT: X-Ray Diffraction Computed Tomography reconstruction tool

PyXRDCT is an X-Ray Diffraction Computed Tomography reconstruction tool.
The library gives functions for reconstructing data from ESRF Bliss h5. So far the library is tested only at ID11 beamline.

## Installation

### Cloning

The source code is available here by running the command:

	git clone https://github.com/poautran/PyXRDCT.git

Then run the following commands to install it:

	cd PyXRDCT
	python3 -m pip install .

The package should be working, to test it run in your python env:

	import PyXRDCT

## Dependencies and OS

The library is currently tested on Python 3.8.10 over Ubuntu 20.04. Further version will be tested in the future.

### Requirements :

	* numpy 		- 	http://www.numpy.org
	* scipy 		- 	http://www.scipy.org
	* matplotlib 		- 	http://matplotlib.sourceforge.net/
	* h5py	    		-  	http://www.h5py.org/
	* fabio			-	https://github.com/silx-kit/fabio
	* pyFAI			-	https://pyfai.readthedocs.io/
	* Others, not detailed yet

## Documentation

The user manual is under construction. For an early start, you can try to run the demo Notebook in PyXRDCT/resources. Functions and other utilities are in PyXRDCT/nmutils/utils. Examples of raw data are located in the resources folder in PyXRDCT/resources.



