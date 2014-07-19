###############################################################################
# Spherical coordinate systems and projections
###############################################################################

from math import *
import numpy as np

###############################################################################
rad = pi/180.0; twopi = 2.0*pi; halfpi = pi/2.0
deg = 180.0/pi; arcmin = deg/60.0; arcsec = arcmin/60.0
###############################################################################

def iau2ang(l,b):
    """Convert IAU galactic coordinates (in degrees) to sperical coordinates"""
    theta = (90.0-b)*rad; phi = l*rad; return (theta,phi)

def ang2iau(theta,phi):
    """Convert sperical coordinates to IAU galactic coordinates (in degrees)"""
    l = phi*deg; b = 90.0 - theta*deg; return (l,b)

def ang2vec(theta,phi):
    """Convert spherical coordinates (in radians) to a cartesian vector"""
    return np.array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])

def iau2vec(l,b):
    """Convert IAU galactic coordinates (in degrees) to a cartesian vector"""
    theta, phi = iau2ang(l,b); return ang2vec(theta,phi)

###############################################################################
from matplotlib.projections.geo import GeoAxes

def latitude_hpx2mpl(x):
    """Wrap latitude from healpix [0,pi] to matplotlib [-pi/2,pi/2] convention"""
    return halfpi - x

def latitude_mpl2hpx(x):
    """Wrap latitude from matplotlib [-pi/2,pi/2] to healpix [0,pi] convention"""
    return halfpi - x

def longitude_hpx2mpl(x):
    """Wrap longitude from healpix [0,2pi] to matplotlib [-pi,pi] convention"""
    return twopi - x if x > pi else (-x if x != 0 else 0)

def longitude_mpl2hpx(x):
    """Wrap longitude from matplotlib [-pi,pi] to healpix [0,2pi] convention"""
    return twopi - x if x > 0  else (-x if x != 0 else 0)

class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts longitude labelling in Mollweide plots from -180,180 to 0-360"""
    def __call__(self, x, pos=None):
        return GeoAxes.ThetaFormatter.__call__(self, longitude_mpl2hpx(x), pos)

def rho(theta):
    return 2.0*np.sin(theta/2.0)

lat = np.vectorize(latitude_hpx2mpl)
lon = np.vectorize(longitude_hpx2mpl)
