###############################################################################
# Spherical coordinate systems and projections
###############################################################################

from math import *
import numpy as np

###############################################################################

def ang2vec(theta,phi):
    """Convert spherical coordinates (in radians) to a cartesian vector"""
    return np.array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])

def iau2vec(l,b):
    """Convert IAU galactic coordiates (in degrees) to a cartesian vector"""
    theta = (90.0-b) * pi/180.0; phi = l * pi/180.0; return ang2vec(theta,phi)
