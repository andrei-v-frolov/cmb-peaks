#!/usr/bin/env python
# calculate min and max averages a la Larson&Wandelt
# usage: minmax.py [peak data]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from math import *
import numpy as np

###############################################################################
# load data and compute averages
###############################################################################

assert len(sys.argv) > 1, "usage: minmax.py [peak data]"

peaks = np.loadtxt(sys.argv[1]); n,m = peaks.shape

assert m > 3, "peak data malformed"

mins = [peaks[i,2] for i in range(n) if peaks[i,3] < 0]; nmin = len(mins)
maxs = [peaks[i,2] for i in range(n) if peaks[i,3] > 0]; nmax = len(maxs)

print -np.sum(mins)/nmin, np.sum(maxs)/nmax, nmin, nmax, n
