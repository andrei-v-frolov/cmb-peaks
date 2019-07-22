#!/usr/bin/env python
# Plot HEALPix maps in Mollweide projection
# usage: plot-map <plot format, e.g. QU/IQU> <data.fits> [output.pdf]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys, warnings
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from math import *
import numpy as np
import healpy as hp
from sys import argv, stdin

# FITS I/O libraries
try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits
import numpy as np


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
assert len(argv) > 2, 'usage: plot-map <plot format, e.g. QU/IQU> <data.fits> [output.pdf]'

plot,maps = argv[1].split('/'); data = argv[2]
file = argv[3] if len(argv) > 3 else None


###############################################################################
# read map data
###############################################################################

# you wish it would work like this...
#smica = hp.read_map(base + 'dx12_v3_smica_cmb_ieb_040a_0256.fits')

def read_map(file, path='.', format='IEB'):
	data = pyfits.open(os.path.join(path, file))[1]
	
	# parse metadata
	nside = data.header['NSIDE']
	order = data.header['ORDERING']
	
	nested = None
	if (order == 'RING'): nested = False
	if (order == 'NESTED'): nested = True
	
	# override column names
	for (i,c) in enumerate(format):
		data.columns[i].name = c
	
	return {c: np.ma.masked_invalid(data.data[c].flatten()) for c in format}, nside, nested

M,nside,nested = read_map(data, format=maps)


###############################################################################
# create the plots
###############################################################################

os.environ['matplotlib.backend'] = 'macosx' if file is None else 'pdf'

# configure matplotlib options
from pltconfig import *
from colormaps import *
import matplotlib

n = len(plot); height = width/(2*n) + 0.3*width
a = min([M[p].min() for p in plot])
b = max([M[p].max() for p in plot])
t = matplotlib.ticker.MaxNLocator().tick_values(a,b); a = t[0]; b = t[-1]

fig = plt.figure(figsize=(cm2inch(width), cm2inch(height)), frameon=fill)

transparent = '#00000000'; grey = '#808080A0'; planck_cmap.set_under(transparent)
c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", [grey, transparent])
c.set_bad('red'); c.set_over('black'); c.set_under(transparent)

for (i,p) in enumerate(plot):
	hp.mollview(M[p].filled(hp.UNSEEN), nest=nested, sub=(1,n,i+1), title="$"+p+"$", cbar=False, min=a, max=b, cmap=planck_cmap, xsize=2400)
	hp.graticule(dmer=360,dpar=360,alpha=0); bbox = plt.gca().get_position(); bbox.x0 -= 0.01; bbox.x1 -= 0.01; plt.gca().set_position(bbox)

###############################################################################
# color bars
###############################################################################

ax = plt.axes([0.1,1.7/height,0.8,0.25/height])
bar = matplotlib.colors.Normalize(vmin=a, vmax=b)
matplotlib.colorbar.ColorbarBase(ax, cmap=planck_cmap, norm=bar, orientation='horizontal', ticks=t).set_label(r'$\mu$K')

###############################################################################

#plt.tight_layout()
plt.show()

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
