#!/usr/bin/env python
# Calculate Kolmogorov-Smirnov deviation from reference distribution
usage = 'ks-test <peaks.dat> <reference.dat> [scaling factor]'

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


###############################################################################
# import peak data
###############################################################################

# parse comand line arguments
assert len(argv) > 2 and len(argv) < 5, usage

# peak data format: theta, phi, value, kind;
# sims data format: value, 9 percentile brackets;

# load peak data, and bracketed sim CDFs
peaks = np.loadtxt(argv[1] if argv[1] != '-' else stdin)
sims  = np.loadtxt(argv[2] if argv[2] != '-' else stdin)

# fudge simulation data if requested
if (len(argv) > 3):
    FUDGE = float(argv[3]); sims[:,0] *= FUDGE

# form peak CDF
x, f, n = makecdf(peaks[:,2])
y = np.vectorize(clint1d(sims[:,0], sims[:,5]))(x)

# Kolmogorov-Smirnov deviation
K = sqrt(n)+0.12+0.11/sqrt(n)

print K*np.amax(np.fabs(f-y))
