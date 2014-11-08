#!/usr/bin/env python
# Calculate significance of the coldest peak (using simulation bootstrap)
# usage: bootstrap <peak data> <coldest peak distribution> [scaling factor]

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
    
    cdf = interp1d(x, f, kind='linear')
    
    return log10(cdf(v))


###############################################################################
# import peak data and do boostrap analysis of significance
###############################################################################
assert len(argv) > 2 and len(argv) < 5, 'usage: bootstrap <peak data> <coldest peak distribution> [scaling factor]'

# parse kernel info
kernel, fwhm = parse(argv[1])

# load peak data and statistics
DATA = np.loadtxt(argv[1]); p = DATA[0,2]
SIMS = np.loadtxt(argv[2]); X = SIMS[:,2]

# scale sims if fudge factor is supplied
if len(argv) == 4: X = X * float(argv[3])

# baseline significance
sign = logutp(p, X); samples = 10000

# bootstrap
###############################################################################
# run the job in parallel, if job control is available
###############################################################################

# parallel worker routine: boostrap significance variance
def varsign(i): return logutp(p, bootstrap(X)) - sign

try:
    from joblib import Parallel, delayed
    from multiprocessing import cpu_count
    
    ds = np.array(Parallel(n_jobs=cpu_count())(delayed(varsign)(i) for i in range(samples)))
except ImportError:
    ds = np.array([varsign(i) for i in range(samples)])

# output significance
print fwhm, sign, sqrt(np.sum(ds*ds)/samples)
