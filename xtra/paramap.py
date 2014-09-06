#!/usr/bin/env python
# 
# usage: paramap [input datafile] [diask diameter in deg, default is 30]

###############################################################################
# import libraries
###############################################################################

from sys import argv, stdin, stdout
from peakstats import *

# parse arguments
file = argv[1] if (len(argv) > 1 and argv[1] != '-') else stdin
fwhm = float(argv[2])*pi/180.0 if (len(argv) > 2) else pi/6.0

###############################################################################
# import peak data
###############################################################################

# peak data format: theta, phi, value, kind (-1 => min, +1 => max);
peaks = np.loadtxt(file)
pixel = np.loadtxt('pixels-nest-64.dat')

npks = peaks.shape[0]
npix = pixel.shape[0]

# fit full-sky peak CDF with Gaussian random peak distribution
fullsky_fit, fullsky_cov = cdf_fit(np.sort(peaks[:,2]), np.linspace(0.0, 1.0, npks))


###############################################################################
# scan peak distribution parameters in a local cap
###############################################################################

def ang2vec(theta, phi):
    """Convert spherical coordinates (in radians) to a cartesian vector"""
    return np.array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])

def cap(theta, phi, fwhm=8.0*pi/180.0):
    """Extract peaks in a disc around a given direction"""
    
    pole = ang2vec(theta,phi)
    qcut = cos(fwhm/2.0)
    data = []
    
    for i in range(npks):
        v = ang2vec(peaks[i,0], peaks[i,1])
        if (np.dot(v, pole) > qcut):
            data.append(peaks[i,2])
    
    return data

def ksdiff(x,f,fit):
    """Calculate Kolmogorov-Smirnov deviation from Gaussian CDF fit"""
    
    # Kolmogorov-Smirnov deviation
    N = len(x); K = sqrt(N)+0.12+0.11/sqrt(N)
    
    # evaluate best fit CDF and fit variance
    gamma, sigma, alpha = fit; y = CDF(x, gamma, sigma, alpha)
    
    return K*max(np.abs(f-y))

def lparams(p):
    """Calculate and print local peak distribution parameters in a disk centered on pixel p"""
    
    theta = pixel[p,0]; phi = pixel[p,1]
    x, f, n = makecdf(cap(theta, phi, fwhm))
    
    print theta, phi, n, ksdiff(x,f,fullsky_fit),
    
    try:
        # fit Gaussian random peak distribution
        fit, cov = cdf_fit(x,f); gamma, sigma, alpha = fit
        print gamma, sigma, alpha, ksdiff(x,f,fit)
    except:
        # not converged, do nothing
        print ""
    
    stdout.flush()


###############################################################################
# run the job in parallel, if job control is available
###############################################################################
try:
    from joblib import Parallel, delayed
    from multiprocessing import cpu_count
    
    Parallel(n_jobs=cpu_count(), verbose=15)(delayed(lparams)(i) for i in range(npix))
except:
    for i in range(npix):
        lparams(i)
