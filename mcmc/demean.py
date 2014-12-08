#!/usr/bin/env python
# Remove mean from peak distribution for further analysis
# usage: demean <peak data>

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from peakstats import *
from sys import argv, stdout
from os.path import basename

###############################################################################
# import peak data and remove distribution mean
###############################################################################
assert len(argv) == 2, 'usage: demean <peak data>'

# load peak data and statistics
DATA = np.loadtxt(argv[1])

# write de-meaned data to stdout
DATA[:,2] = demean(DATA[:,2])
np.savetxt(stdout, DATA, fmt=['%24.16g','%24.16g','%24.16g','%i'])
