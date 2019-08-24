#!/usr/bin/env python
# Plot 3D vector field spectra
# usage: plot-vspectra <vector.fits> [output.pdf]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys, warnings
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from mapio import *
import numpy as np
from sys import argv, stdin


###############################################################################
# plot style options
###############################################################################

width =  8.0	# figure width, in centimeters
aspect = 3/4.	# figure aspect ratio (default is 4:3)
fill = False	# produce transparent plots if false
grid = True	# do we want to render the plot grid?
legend = True	# do we want to include the legend?


###############################################################################
# parse arguments
###############################################################################
assert len(argv) > 1, 'usage: plot-vspectra <vector.fits> [output.pdf]'

file = argv[2] if len(argv) > 2 else None

###############################################################################
# create the plots
###############################################################################

os.environ['matplotlib.backend'] = 'macosx' if file is None else 'pdf'

# configure matplotlib options
from pltconfig import *

fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)


###############################################################################

M = hp.read_map(argv[1], field=None); lmax = 20
P = np.zeros((3,3,lmax+1)); C = np.zeros((3,lmax+1))

for i in range(3):
	for j in range(3):
		P[i,j,:] = hp.anafast(M[i], M[j], lmax=lmax, pol=False)

for l in range(lmax+1):
	C[:,l] = np.linalg.eigvalsh(P[:,:,l])

l = np.arange(lmax+1)

plt.loglog(l, (C[0,:]+C[1,:]+C[2,:]) * l*(l+1)/(2.0*np.pi), 'k-', label='')

plt.loglog(l, C[0,:] * l*(l+1)/(2.0*np.pi), 'r-', label='', alpha=0.2)
plt.loglog(l, C[1,:] * l*(l+1)/(2.0*np.pi), 'r-', label='', alpha=0.2)
plt.loglog(l, C[2,:] * l*(l+1)/(2.0*np.pi), 'r-', label='', alpha=0.2)

plt.ylabel(r'$\ell(\ell+1)\ C_\ell^{XX}/2\pi$'); plt.xlabel(r'$\ell$')

plt.xlim(1,lmax); plt.ylim(1.0e-6,1.0)

###############################################################################

plt.tight_layout()
plt.show()

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
