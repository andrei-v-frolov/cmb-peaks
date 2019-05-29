#!/usr/bin/env python
# Plot 2D histogram counts
# usage: plot-counts <data.fits> [output.pdf] [downsample]

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
from sys import argv, stdin

# FITS I/O libraries
try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits


###############################################################################
# plot style options
###############################################################################

width = 10.0	# figure width, in centimeters
aspect = 5/4.	# figure aspect ratio (default is 4:3)
fill = False	# produce transparent plots if false
grid = True	# do we want to render the plot grid?
legend = False	# do we want to include the legend?


###############################################################################
# parse arguments
###############################################################################
assert len(argv) > 1, 'usage: plot-counts <data.fits> [output.pdf] [downsample]'

hdu = pyfits.open(argv[1])
file = argv[2] if len(argv) > 2 else None
scale = int(argv[3]) if len(argv) > 3 else 1


###############################################################################
# process the data
###############################################################################

from scipy.signal import medfilt

data = hdu[0].data
ny,nx = data.shape

# antialias using median filter
osy = scale; osx = scale

if (osy != 1 or osx != 1):
	aakernel = [(osy & (~1)) + 1, (osx & (~1)) + 1]
	data = medfilt(data, aakernel)[::osy,::osx]; ny,nx = data.shape
	print("Antialiased with median kernel %s; output size is (%i,%i) pixels" % (aakernel,nx,ny))

# data extent and meshes
x0 = hdu[0].header['CRVAL1']
y0 = hdu[0].header['CRVAL2']
dx = hdu[0].header['CDELT1']
dy = hdu[0].header['CDELT2']

x = x0 + np.arange(nx)*dx*osx
y = y0 + np.arange(ny)*dy*osy

X, Y = np.meshgrid(x, y)

hdu.close()

###############################################################################
# create the plots
###############################################################################

os.environ['matplotlib.backend'] = 'macosx' if file is None else 'pdf'

# configure matplotlib options
from pltconfig import *

import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)

# plot image data
#c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", ["blue", "white", "black"])
img = plt.imshow(data, extent=(x[0],x[-1],y[0],y[-1]), origin='lower', aspect=(osx*dx)/(osy*dy), cmap=matplotlib.cm.Oranges, interpolation='none')

# plot contours
#plt.contour(X, Y, data, 32, cmap=matplotlib.cm.jet)

if grid: plt.grid(True, which="both", axis="both", linestyle='-', linewidth=0.2)
if legend: plt.legend(loc='upper right', ncol=2, frameon=fill or grid)

# make colorbar match the plot
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=-0.4)
plt.colorbar(img, cax=cax)

###############################################################################

plt.tight_layout()
plt.show()

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
