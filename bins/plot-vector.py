#!/usr/bin/env python
# Plot 3D vector field on a shaded sphere
# usage: plot-vector <basemap.fits> <vectors.fits> [output.pdf]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys, warnings
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from mapio import *
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
assert len(argv) > 2, 'usage: plot-vector <basemap.fits> <vectors.fits> [output.pdf]'

file = argv[3] if len(argv) > 3 else None

###############################################################################
# create the plots
###############################################################################

os.environ['matplotlib.backend'] = 'macosx' if file is None else 'pdf'

# configure matplotlib options
from pltconfig import *
from colormaps import *
import matplotlib

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(cm2inch(width), cm2inch(width)), frameon=fill)
ax = fig.gca(projection='3d'); ax.set_aspect("equal"); vside = 4; grid = 128

###############################################################################

M = hp.read_map(argv[1]); nside = hp.get_nside(M)

theta = np.linspace(np.pi, 0, grid)
phi = np.linspace(-np.pi, np.pi, grid)
phi,theta = np.meshgrid(phi, theta)
gridmap = M[hp.ang2pix(nside, theta, phi)]
colors = planck_cmap(gridmap/gridmap.max())

x = np.sin(theta)*np.cos(phi)
y = np.sin(theta)*np.sin(phi)
z = np.cos(theta)

# plot the surface colored with the base map
ax.plot_surface(x, y, z, facecolors=colors, shade=False, rcount=grid, ccount=grid)

###############################################################################

M,nside,nested = read_map(argv[2], format='xyz')

x,y,z = hp.pix2vec(vside, [range(12*vside**2)], nest=nested)

u = hp.ud_grade(M['x'], vside, order_in='NESTED' if nested else 'RING')
v = hp.ud_grade(M['y'], vside, order_in='NESTED' if nested else 'RING')
w = hp.ud_grade(M['z'], vside, order_in='NESTED' if nested else 'RING')

c = np.sqrt(u*u+v*v+w*w).flatten(); c = c/c.max()
c = np.concatenate((c, np.repeat(c, 2)))
c = plt.get_cmap('winter')(c)

# plot vectors pointing into sphere colored with field strength
ax.quiver(x, y, z, u, v, w, colors=c, pivot='tip', length=0.2, normalize=True, alpha=0.5)

###############################################################################

plt.tight_layout()
plt.show()

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
