#!/usr/bin/env python
# Plot 1-point PDF of HEALPix maps, outputting L-moments for diagnostics
# usage: plot-1pt-pdf <plot format, e.g. QU/IQU/uK> <data.fits> [output.pdf]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys, warnings
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from mapio import *
from peakstats import *
from sys import argv, stdin


###############################################################################
# plot style options
###############################################################################

width = 16.0	# figure width, in centimeters
fill = False	# produce transparent plots if false
grid = True	# do we want to render the plot grid?
legend = True	# do we want to include the legend?

###############################################################################
# parse arguments
###############################################################################
assert len(argv) > 2, 'usage: plot-1pt-pdf <plot format, e.g. QU/IQU/uK> <data.fits> [output.pdf]'

plot,maps,units = (argv[1]+'/').split('/')[:3]
M,nside,nested = read_map(argv[2], format=maps)
file = argv[3] if len(argv) > 3 else None

autoconvert_polarization(M, want=plot)

def makepdf(x, F):
    """Make PDF P(x) of supplied distribution data by finite difference"""
    
    n = len(x); q = [(x[i]+x[i-1])/2 for i in range(1,n)]
    P = [(F[i]-F[i-1])/(x[i]-x[i-1]) for i in range(1,n)]
    
    return q, P

###############################################################################
# create the plots
###############################################################################

os.environ['matplotlib.backend'] = 'macosx' if file is None else 'pdf'

# configure matplotlib options
from pltconfig import *
from colormaps import *
import matplotlib

n = len(plot); height = 0.55*width
a = min([M[p].min() for p in plot])
b = max([M[p].max() for p in plot])
t = matplotlib.ticker.MaxNLocator().tick_values(a,b); a = t[0]; b = t[-1]
X = np.linspace(a, b, 101); A = np.inf; B = 0.0

fig = plt.figure(figsize=(cm2inch(width), cm2inch(height)), frameon=fill)

for (i,p) in enumerate(plot):
	x, F, n = makecdf(M[p].compressed()); L1, L2, L3, L4 = lmoments(x, F)
	print argv[2], p, L2, L3/L2, L4/L2 - 0.12260171954089094743716661166353633
	
	# Gaussian distribution with mean = L1 and sigma = L2/sqrt(pi)
	G = np.exp(-(X-L1)**2/L2**2/(2*np.pi))/np.sqrt(2.0*np.pi**2)/L2; plt.semilogy(X, G,'--')
	
	# decimated PDF with bins selected uniformly in CDF arc length
	x,y = decimate(x, F, 51); x,y = makepdf(x, y); plt.semilogy(x, y, 'o', label="$"+p+"$", ms=3.0)
	
	# sampled PDF limits
	A = min(min(y),A); B = max(max(y),B)

if grid: plt.grid(True, which="both", axis="both", linestyle='-', linewidth=0.2)
if legend: plt.legend(loc='upper left', ncol=2, frameon=fill or grid)

plt.xlim(a,b); plt.xticks(t); plt.xlabel(units); plt.ylim(A/2,B*2)

###############################################################################

plt.tight_layout()
plt.show()

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
