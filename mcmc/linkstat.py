#!/usr/bin/env python
# build peak merger tree
# usage: links.py [peak data files, ordered from fine to coarse scales]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from math import *
import numpy as np
import networkx as nx
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
    #print i, argv[-i]
    
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

#print 'Total %i nodes, %i edges' % (len(G.nodes()), len(G.edges()))

###############################################################################
# tabulate and output degree histogram of a merger tree graph
###############################################################################

bins = 16; count = np.zeros(bins, np.int)

for i in G.degree().values():
    count[i if i < bins else -1] += 1

for i in count:
    print i,
