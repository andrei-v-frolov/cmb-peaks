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

def read_map(file, path='.', hdu=1, format='IQU'):
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

###############################################################################
# convert polarization formats on the fly
###############################################################################

def QU2EB(map, QU='QU', EB='EB', nested=False, mask=None):
	Q = map[QU[0]]; U = map[QU[1]]; I = np.zeros_like(Q)
	if nested: Q = hp.pixelfunc.reorder(Q, n2r=True)
	if nested: U = hp.pixelfunc.reorder(U, n2r=True)
	nside = hp.pixelfunc.get_nside(Q)
	_,ae,ab = hp.sphtfunc.map2alm([I,Q,U], pol=True)
	E = hp.sphtfunc.alm2map(ae, nside, pol=False)
	B = hp.sphtfunc.alm2map(ab, nside, pol=False)
	if nested: E = hp.pixelfunc.reorder(E, r2n=True)
	if nested: B = hp.pixelfunc.reorder(B, r2n=True)
	map[EB[0]] = np.ma.array(E, mask=mask)
	map[EB[1]] = np.ma.array(B, mask=mask)

def EB2QU(map, EB='EB', QU='QU', nested=False, mask=None):
	E = map[EB[0]]; B = map[EB[1]]
	if nested: E = hp.pixelfunc.reorder(E, n2r=True)
	if nested: B = hp.pixelfunc.reorder(B, n2r=True)
	ae,ab = hp.sphtfunc.map2alm([E,B], pol=False)
	ai = np.zeros_like(ae); nside = hp.pixelfunc.get_nside(E)
	_,Q,U = hp.sphtfunc.alm2map([ai,ae,ab], nside, pol=True)
	if nested: Q = hp.pixelfunc.reorder(Q, r2n=True)
	if nested: U = hp.pixelfunc.reorder(U, r2n=True)
	map[QU[0]] = np.ma.array(Q, mask=mask)
	map[QU[1]] = np.ma.array(U, mask=mask)

def autoconvert_polarization(map, want=None, nested=False, mask=None):
	for c in want:
		if ((c in 'EB') and not(c in map) and ('Q' in map) and ('U' in map)):
			print "Converting QU->EB"; QU2EB(map, 'QU', 'EB', nested, mask)
		if ((c in 'QU') and not(c in map) and ('E' in map) and ('B' in map)):
			print "Converting EB->QU"; EB2QU(map, 'EB', 'QU', nested, mask)
		if ((c in 'eb') and not(c in map) and ('q' in map) and ('u' in map)):
			print "Converting qu->eb"; QU2EB(map, 'qu', 'eb', nested, mask)
		if ((c in 'qu') and not(c in map) and ('e' in map) and ('b' in map)):
			print "Converting eb->qu"; EB2QU(map, 'eb', 'qu', nested, mask)
