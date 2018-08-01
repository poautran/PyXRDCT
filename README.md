# PyXRDCT: X-Ray Diffraction Computed Tomography reconstruction functions

PyXRDCT provides X-Ray Diffraction Computed Tomography reconstruction functions.
The library gives functions for correcting diffraction tomography sinograms based on simple mathematics such as center of mass or outliers values.

## Installation

### With PIP

The library is currently under TestPyPI at this link:

	https://test.pypi.org/project/PyXRDCT/

To install the library you can run the command:

	pip install --user --index-url https://test.pypi.org/simple/ PyXRDCT

### From source code

To install PyXRDCT you can download the source code at this url:

	https://github.com/poautran/PyXRDCT/archive/master.zip

Then run the following command to unzip the file:

	unzip PyXRDCT-master.zip

The files and in PyXRDCT-master, then install the library with:

	cd PyXRDCT-master
	python setup.py sdist bdist_wheel

## Dependencies and OS

The library is currently tested on Python 2.7 over Ubuntu 16.04. Further version will be tested in the future.

### Requirements :

	* numpy 		- 	http://www.numpy.org
	* scipy 		- 	http://www.scipy.org
	* matplotlib 		- 	http://matplotlib.sourceforge.net/
	* h5py	    		-  	http://www.h5py.org/

## Documentation

The user manual is under construction. For an early start, main code to be run and modified is in PyXRDCT/PyXRDCT. Functions and other utilities are in PyXRDCT/nmutils/utils. Examples of raw data are located in the resources folder in PyXRDCT/nmutils/resources.



