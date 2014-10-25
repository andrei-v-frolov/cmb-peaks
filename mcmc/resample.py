#!/usr/bin/env python
# Calculate percentile brackets of a CDF using simulated samples
# usage: resample <reference> [MC distribution samples]

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

###############################################################################

# reference grid
def reference(X, xlim=[-6,6], pts=64):
    """Create reference sampling grid covering data set X"""
    
    # extract peak CDF
    x, f, n = makecdf(X)
    fit, cov = cdf_fit(x,f)
    gamma, sigma, alpha = fit
    
    return sigma*np.linspace(xlim[0], xlim[1], pts)

# resample CDF
def resample(X, Y):
    """Resample CDF specified by a distribution X onto points Y"""
    
    x, f, n = makecdf(X); cdf = interp1d(x, f, kind='linear')
    clipcdf = lambda v: 0.0 if v < x[0] else (1.0 if v > x[-1] else cdf(v))
    
    return np.vectorize(clipcdf)(Y)

# percentile brackets
def brackets(X, Y):
    """Evaluate percentile brackets Y of a distribution X"""
    
    x, f, n = makecdf(X); cdf = interp1d(f, x, kind='linear')
    
    return np.vectorize(cdf)(Y)


###############################################################################
# import peak data and marginalize CDFs
###############################################################################

assert len(argv) > 2, 'Not enough arguments, at least reference should be supplied'
assert len(argv) > 4, 'Not enough simulations, only %i supplied' % (len(argv)-2)

# load reference point data
DATA = np.loadtxt(argv[1])
x = reference(DATA[:,2])
n = len(x); nsims = len(argv)-2

# effective sigma brackets to output
sigmas = np.linspace(-3.0,3.0,7)
sigmas = np.vectorize(lambda nu: (1.0+erf(nu/sqrt(2.0)))/2.0)(sigmas)
sigmas = np.concatenate(([0.0],sigmas,[1.0]))


###############################################################################
# run the job in parallel, if job control is available
###############################################################################

# parallel worker routines: resample to reference grid and CDF brackets
def readsim(i): return resample(np.loadtxt(argv[i+2])[:,2], x)
def dumpcdf(i): return brackets(SIMS[:,i], sigmas)

try:
    from joblib import Parallel, delayed
    from multiprocessing import cpu_count
    
    SIMS = np.array(Parallel(n_jobs=cpu_count())(delayed(readsim)(i) for i in range(nsims)))
    BRKS = np.array(Parallel(n_jobs=cpu_count())(delayed(dumpcdf)(i) for i in range(n)))
except ImportError:
    SIMS = np.array([readsim(i) for i in range(nsims)])
    BRKS = np.array([dumpcdf(i) for i in range(n)])

# output CDF brackets on resampled grid
for i in range(n):
    print x[i], ("%.16g " * len(sigmas)) % tuple(BRKS[i,:])
