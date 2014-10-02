#!/usr/bin/env python
# Calculate significance of coldest peak 
# usage: coldspots [peak data files]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from peakstats import *
from sys import argv, stdin
from os.path import basename

# filename handling
def parse(path):
    """Parse peak data file name, extracting kernel type and FWHM"""
    
    name, dot, ext = basename(path).partition('.')
    lmom, kernel, radius = name.split('-')
    
    return kernel, 2.0*float(radius)

###############################################################################
# import peak data and build merger tree
###############################################################################

for i in range(1,len(argv)):
    # parse kernel info
    kernel, fwhm = parse(argv[i])
    
    # read peak data: theta, phi, value, kind (-1 => min, +1 => max);
    peaks = np.loadtxt(argv[i]); n, m = peaks.shape
    assert m > 3, 'not enough peak data columns in ' + argv[i]
    
    # form peak CDF
    x, f, n = makecdf(peaks[:,2])
    
    # try fitting local Gaussian random peak distribution
    try:
        fit, cov = cdf_fit(x,f)
        
        coldsign = marginalize(lambda p: log10(1.0 - (1.0-CDF(x[0], p[0], p[1], p[2]))**n), fit, cov)
        hotsign = marginalize(lambda p: log10(1.0 - CDF(x[-1], p[0], p[1], p[2])**n), fit, cov)
        
        print fwhm, coldsign[0], coldsign[1]
    except Exception:
        # not converged, do nothing
        pass
