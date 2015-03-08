#!/usr/bin/env python
# Calculate significance of the coldest peak (using simulation bootstrap)
usage = 'bootstrap [log]{ltp|utp} <data[:column]> <distribution[:column]> [scaling factor]'

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


###############################################################################
# callable significance routines, all have exactly two arguments
###############################################################################

def ltp(v, X):
    """Lower tail probability to find value x in distribution X (i.e. #X < v)"""
    x, f, n = makecdf(X); return lut1d(v, x, f)

def utp(v, X):
    """Upper tail probability to find value x in distribution X (i.e. #X > v)"""
    x, f, n = makecdf(X); return lut1d(v, x, 1.0-f)

def logltp(v, X):
    """Log-lower tail probability to find value x in distribution X (i.e. #X < v)"""
    x, f, n = makecdf(X); f[0] = 1.0/n; return lut1d(v, x, np.log10(f))

def logutp(v, X):
    """Log-upper tail probability to find value x in distribution X (i.e. #X > v)"""
    x, f, n = makecdf(X); f[-1] = (n-1.0)/n; return lut1d(v, x, np.log10(1.0-f))


###############################################################################
# import peak data and do boostrap analysis of significance
###############################################################################

# parse comand line arguments
assert len(argv) > 3 and len(argv) < 6, usage
assert argv[1] in locals(), usage

significance = locals()[argv[1]]
data,dcol = (argv[2]+':2').split(':')[0:2]
sims,scol = (argv[3]+':2').split(':')[0:2]
FUDGE = float(argv[4]) if len(argv) > 4 else None

# parse kernel info
try:
    kernel, fwhm = parse(data)
except ValueError:
    kernel, fwhm = parse(sims)

# load peak data and statistics
DATA = np.loadtxt(data if data != '-' else stdin); p = DATA[0,int(dcol)] if DATA.ndim > 1 else (DATA[int(dcol)] if DATA.ndim > 0 else DATA)
SIMS = np.loadtxt(sims if sims != '-' else stdin); X = SIMS[:,int(scol)] if SIMS.ndim > 1 else SIMS

# scale sims if fudge factor is supplied
if FUDGE: X *= FUDGE

# baseline significance
sign = significance(p, X); samples = 10000


###############################################################################
# bootstrap, run the job in parallel if job control is available
###############################################################################

# parallel worker routine: boostrap significance variance
def varsign(i): return significance(p, bootstrap(X)) - sign

try:
    from joblib import Parallel, delayed
    from multiprocessing import cpu_count
    
    ds = np.array(Parallel(n_jobs=cpu_count())(delayed(varsign)(i) for i in range(samples)))
except ImportError:
    ds = np.array([varsign(i) for i in range(samples)])

# output significance
print fwhm, sign, sqrt(np.sum(ds*ds)/samples)
