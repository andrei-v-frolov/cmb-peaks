#!/usr/bin/env python
# build peak merger tree
# usage: links.py [peak data files, ordered from fine to coarse scales]

###############################################################################
# import libraries
###############################################################################

from math import *
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt
from matplotlib import gridspec

from colormaps import *
from spherical import *

###############################################################################
# file handling
###############################################################################

from sys import argv, stdin
from os.path import basename

def parse(path):
    """Parse peak data file name, extracting kernel type and FWHM"""
    
    name, dot, ext = basename(path).partition('.')
    lmom, kernel, radius = name.split('-')
    
    return kernel, 2.0*float(radius)

###############################################################################
# import peak data and build merger tree
###############################################################################

coarse = {}
scales = []

G = nx.Graph()

def closest(x, data):
    """Find a node closest to x in a scale space layer"""
    
    # unpack scale space data
    kernel, fwhm, nodes = data
    
    # evaluate distance measure
    dx = nodes - x
    ds = (dx*dx).sum(axis=1)
    idx = np.argmin(ds)
    
    # return closest node tuple
    return (kernel,fwhm,idx) if ds[idx] < (fwhm/arcmin)**2/4.0 else None

for i in range(1,len(argv)):
    print i, argv[-i]
    
    # parse kernel info
    kernel, fwhm = parse(argv[-i]); scales.append(fwhm)
    
    # read peak data: theta, phi, value, kind (-1 => min, +1 => max);
    peaks = np.loadtxt(argv[-i]); n, m = peaks.shape
    assert m > 3, 'not enough peak data columns in ' + argv[-i]
    
    # node vector space
    nodes = np.empty([n,5])
    
    # index peaks according to their value
    idx = np.argsort(peaks[:,2],  kind='mergesort')
    
    # build merger tree
    for rank in range(n):
        p = peaks[idx[rank],0:2]
        v = peaks[idx[rank],2]
        k = peaks[idx[rank],3]
        f = rank/(n-1.0)
        
        # create a node vector
        node = np.append(ang2vec(p[0], p[1]), [f, k])
        
        # various projections
        moll = [lon(p[1]), lat(p[0])]
        lamn = [p[1], rho(p[0])]
        lams = [p[1], rho(pi-p[0])]
        tree = [f, fwhm]
        
        # node size
        size = 40.0*(fwhm/800.0)**2
        
        # add current node to the graph
        G.add_node((kernel,fwhm,rank),
            pos=p, val=v, rank=f, kind=k, size=size,
            mollweide=moll, north=lamn, south=lams, tree=tree)
        
        # link to the closest coarse node
        if (kernel in coarse):
            link = closest(node, coarse[kernel])
            
            if (link):
                G.add_edge((kernel,fwhm,rank), link)
        
        # stash node vector
        nodes[rank,:] = node
    
    # done; update coarse layer
    coarse[kernel] = [kernel, fwhm, nodes]

print 'Total %i nodes, %i edges' % (len(G.nodes()), len(G.edges()))

###############################################################################
# spherical projections
###############################################################################

def mollweide_setup(sky, grid=True, latgrid=30, longrid=60, latmax=90):
    """Setup Mollweide projection plot"""
    
    # graticule
    sky.grid(grid)
    sky.set_latitude_grid(latgrid)
    sky.set_longitude_grid(longrid)
    sky.set_longitude_grid_ends(latmax)
    
    #sky.xaxis.set_ticklabels([])
    sky.yaxis.set_ticklabels([])
    sky.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))

def lambert_setup(sky, origin='S', direction=-1, grid=True, latgrid=30, longrid=60, latmax=90):
    """Setup Lambert projection plot"""
    
    tgrid = np.arange(0,360,longrid)
    rgrid = rho(np.arange(latmax,0,-latgrid)*rad)
    
    sky.grid(grid)
    sky.set_theta_zero_location(origin)
    sky.set_theta_direction(direction)
    sky.set_thetagrids(tgrid)
    sky.set_rgrids(rgrid, [])
    sky.set_rmax(rgrid[0])

###############################################################################
# draw merger tree
###############################################################################
fig = plt.figure(figsize=(16,8), frameon=False)

layout = gridspec.GridSpec(2, 3, height_ratios=[1, 1], width_ratios=[1, 2, 1],
               left=0.03, right=0.97, bottom=0.03, top=0.97, wspace=0.0, hspace=0.15)

# views of the sky
sky = plt.subplot(layout[0,1], projection='mollweide')
north = plt.subplot(layout[0,0], projection='polar')
south = plt.subplot(layout[0,2], projection='polar')

# merger tree
tree = plt.subplot(layout[1,:], yscale='log')
plt.xlim([-0.01,1.01]); plt.ylim([min(scales)/1.05,max(scales)*1.05])

# draw graph
rank = [G.node[i]['rank'] for i in nx.nodes_iter(G)]
size = [G.node[i]['size'] for i in nx.nodes_iter(G)]

nx.draw_networkx_nodes(G, nx.get_node_attributes(G,'mollweide'), ax=sky,
    with_labels=False, node_color=rank, node_size=size, linewidths=0.0, alpha=0.3)
nx.draw_networkx_nodes(G, nx.get_node_attributes(G,'north'), ax=north,
    with_labels=False, node_color=rank, node_size=size, linewidths=0.0, alpha=0.3)
nx.draw_networkx_nodes(G, nx.get_node_attributes(G,'south'), ax=south,
    with_labels=False, node_color=rank, node_size=size, linewidths=0.0, alpha=0.3)
nx.draw_networkx(G, nx.get_node_attributes(G,'tree'), ax=tree,
    with_labels=False, node_color=rank, node_size=size, linewidths=0.0, alpha=0.3)

# finalize setup
mollweide_setup(sky)
lambert_setup(north, 'S', -1)
lambert_setup(south, 'N', +1)

# draw the figure
#plt.savefig('links.png', bbox_inches='tight', pad_inches=0.02, transparent=True)
plt.show()
