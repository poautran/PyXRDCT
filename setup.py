#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyXRDCT",
    version="0.2.0",
    author="Pierre-Olivier Autran",
    author_email="pierre-olivier.autran@esrf.fr",
    description="This packages does provide X-Ray Diffraction Computed Tomography reconstruction functions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/poautran/PyXRDCT",
    packages=setuptools.find_packages(),
    classifiers=(
    	"Development Status :: 1 - Planning",
        "Natural Language :: English",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ),
)
