###############################################################################
# import libraries
###############################################################################

# OS interface
import os

# import libraries
from math import *
import numpy as np
import healpy as hp

# FITS I/O libraries
try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits
import numpy as np


###############################################################################
# read map data
###############################################################################

# you wish it would work like this...
#smica = hp.read_map(base + 'dx12_v3_smica_cmb_ieb_040a_0256.fits')

def read_map(file, path='.', hdu=1, format='IEB'):
	data = pyfits.open(os.path.join(path, file))[hdu]
	
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
