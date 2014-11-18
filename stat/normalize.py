#!/usr/bin/env python
# normalize peak data to sigma units (using Gaussian sky assumption)
# usage: normalize <peak data[:sims[:scale]]>

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
# import peak data
###############################################################################

assert len(argv) > 1, "usage: normalize <peak data[:sims[:scale]]>"

# load peak data, and bracketed sim CDFs
data = argv[1].split(':'); args = len(data)

# peak data format: theta, phi, value, kind;
# sims data format: value, 9 percentile brackets;

if (args == 1):
    peaks = np.loadtxt(data[0] if data[0] != '-' else stdin)
    sims  = None
elif (args <= 3):
    peaks = np.loadtxt(data[0] if data[0] != '-' else stdin)
    sims  = np.loadtxt(data[1] if data[1] != '-' else stdin)
    
    # scale simulation data if fudge factor is supplied
    if (args == 3): sims[:,0] = sims[:,0] * float(data[2])
else:
    raise SystemExit("Don't know what to do with extra arguments to plot data, aborting...")

# parse kernel type and FWHM
kernel,fwhm = parse(data[0])

# form peak CDF
x, f, n = makecdf(peaks[:,2])

try:
    # fit Gaussian random peak distribution
    if (sims is None):
        fit, cov = cdf_fit(x,f)
    else:
        fit, cov = cdf_fit(sims[:,0],sims[:,5])
    gamma, sigma, alpha = fit
    
    # dump normalized peaks
    for i in range(n):
        print peaks[i,0], peaks[i,1], peaks[i,2]/sigma, fwhm
except Exception:
    # not converged, do nothing
    pass
