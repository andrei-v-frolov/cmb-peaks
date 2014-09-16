#!/usr/bin/env python
# Calculate significance of coldest peak 
# usage: coldspots [peak data files]

###############################################################################
# import libraries
###############################################################################

from sys import argv, stdin
from peakstats import *

for i in range(1,len(argv)):
    #print i, argv[i]
    
    # peak data format: theta, phi, value, kind (-1 => min, +1 => max);
    peaks = np.loadtxt(argv[i])
    
    # form peak CDF
    x, f, n = makecdf(peaks[:,2])
    
    # fit Gaussian random peak distribution
    fit, cov = cdf_fit(x,f)
    
    coldsign = marginalize(lambda p: log(1.0 - (1.0-CDF(x[0], p[0], p[1], p[2]))**n), fit, cov)
    hotsign = marginalize(lambda p: log(1.0 - CDF(x[-1], p[0], p[1], p[2])**n), fit, cov)
    
    print argv[i], coldsign[0], coldsign[1]
