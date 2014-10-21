#!/usr/bin/env python
# Calculate significance of the coldest peak (using simulation bootstrap)
# usage: bootstrap [peak data] [coldest peak distribution]

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

# bootstrap
def bootstrap(X):
    """Resample an array using bootstrap method"""
    
    n = len(X); resample = np.floor(np.random.rand(n)*n).astype(int)
    
    return X[resample]

# log-likelihood
def logutp(v, X):
    """Log-upper tail probability to find value x in distribution X"""
    
    x, f, n = makecdf(X)
    
    if (v < x[ 0]): return log10(1.0/len(x))
    if (v > x[-1]): return 0.0
    
    cdf = interp1d(x, f, kind='linear', assume_sorted=True)
    
    return log10(cdf(v))


###############################################################################
# import peak data and do boostrap analysis of significance
###############################################################################

# parse kernel info
kernel, fwhm = parse(argv[1])

# load peak data and statistics
DATA = np.loadtxt(argv[1]); p = DATA[0,2]
SIMS = np.loadtxt(argv[2]); X = SIMS[:,2]

# baseline significance
sign = logutp(p, X); var  = 0.0; samples = 10000

# bootstrap
for i in range(samples):
    ds = logutp(p, bootstrap(X)) - sign; var += ds*ds

# output significance
print fwhm, sign, sqrt(var/samples)
