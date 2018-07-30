#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Importing couple of packages here
import numpy as np 
import sys, os, glob
import matplotlib.pyplot as plt
import h5py

#Definig all the functions required to display results


def createCircularMask(h, w, center=None, radius=None):
    """Creates a circular mask on the picture with b and w as sizes of the circular shape. 
    Radius and center to be defined."""
    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def progression(string,top_value,iterator):
    """Gives an idea of the percentage of progression of the loop. Need to put the string you 
    want before the %, a top value and the name of the loop iterator."""
    sys.stdout.write("\r\x1b[K"+ string +str( '%.2f'%float((float(iterator)+1)*100/top_value))+ "%")
    sys.stdout.flush()
