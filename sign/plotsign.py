#!/usr/bin/env python
# Plot coldest spot significance as a function of filter FWHM
# usage: plotsign [width list (in cm)]

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

# plot style options
aspect = 16/6.  # figure aspect ratio (default is 4:3)
fill = False    # produce transparent plots if false
grid = False    # do we want to render the plot grid?

widths = map(float, sys.argv[1].split(',')) if (len(argv) > 1) else [18.0, 12.0, 8.8]

###############################################################################
# import cold spot data
###############################################################################

try:
    dataset = open(os.path.join(here, '../DATASET'), 'r').readline().strip()
except StandardError:
    dataset = ''

plot = {
    'SSG21': r'\texttt{SSG21} filter' + (', ' + dataset if dataset else ''),
    'SSG42': r'\texttt{SSG42} filter' + (', ' + dataset if dataset else ''),
    'SSG84': r'\texttt{SSG84} filter' + (', ' + dataset if dataset else '')
}

###############################################################################
# create the plots
###############################################################################

# Configure Matplotlib options
from pltconfig import *
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

def plot_setup(ax, xlabel='', ylabel='', xlim=[2,20], ylim=[-4,0]):
    """Setup common plot style parameters"""
    
    # labels
    if (xlabel):
        plt.xlabel(xlabel)
    if (ylabel):
        plt.ylabel(ylabel)
    
    # distance of axis label to tick labels
    ax.yaxis.labelpad = 10*width/17.
    ax.xaxis.labelpad = 10*width/17.
    
    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # grid
    if (grid):
        plt.grid(True, which="major", axis="both")
    
    # axes limits
    plt.ylim(ylim)
    plt.xlim(xlim)
    
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")

# Create the plots
for width in widths:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)
    
    ax = fig.add_subplot(111)
    
    plt.title("Significance of the coldest spot on the sky", y=1.2)
    plot_setup(ax, xlabel='filter kernel size [degrees full width half max]',
        ylabel='probability of Gaussian realization', ylim=[1.0e-4,1.0e0])
    
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    
    xs = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,32]
    plt.xticks(xs, xs)
    
    # effective confidence nu is defined by P(excursion > nu) = P
    plot_setup(ax.twinx(), ylabel=r'effective confidence level [$\sigma$]')
    
    sigmas = np.linspace(0.5,3.5,7)
    ys = np.vectorize(lambda x: log10(erfc(x/sqrt(2.0))/2.0))(sigmas)
    plt.yticks(ys, sigmas)
    
    plt.hlines(ys, 2, 20, linewidth=0.5, linestyle=':', zorder=0)
    plt.fill_between([2,20], ys[5], -4.0, edgecolor='none', color='gold', alpha=0.2, zorder=-5)
    
    for f, l in plot.items():
        try:
            data = np.loadtxt(f) # significance data format: fwhm, logsign, dlogsign;
            plt.errorbar(data[:,0]/60.0, data[:,1], yerr=data[:,2], fmt='-o', alpha=0.8, label=l)
        except StandardError:
            MissingData = f + ": Data file not found!"
            warnings.warn(MissingData)
    
    # legend setup
    if (fill or True):
        leg = plt.legend(loc='lower left', numpoints=1, frameon=True)
        
        if (leg): # remove box around legend
            leg.get_frame().set_edgecolor('none')
            leg.get_frame().set_alpha(.8)
    else:
        plt.legend(loc='lower left', numpoints=1, frameon=False)
    
    # override auto-scale
    plt.ylim([-4,0])
    
    # effective bandwidth axis
    ax.twiny()
    
    BW = 3.97 # for GAUSS filters
    BW = 4.17 # for SSGxx filters
    
    plt.xlabel(r'effective filter bandwidth [$\ell$ at which $b_\ell = b_{max}/2$]')
    ls = [7,8,9,10,11,12,13,15,17,20,25,30,40,50,60,70,80,90,100,110]
    xs = (BW - np.vectorize(log10)(np.array(ls)) - log10(120.0))/log10(1200.0/120.0)
    
    plt.xticks(xs, ls); plt.xlim([0,1])
    
    # reduce white space around figure
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0)
    
    # save to pdf with right bounding box
    base = "sign"
    pdf = "%s.%dmm.pdf" % (base,int(width*10)) if len(widths) > 1 else base+".pdf"
    plt.savefig(pdf, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
