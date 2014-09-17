#!/usr/bin/env python

# configure import path
import os, sys
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# Configure Matplotlib options
from pltconfig import *
from matplotlib.ticker import MaxNLocator

# plot style options
aspect = 3/2.   # figure aspect ratio (default is 4:3)
fill = False    # produce transparent plots if false
grid = False    # do we want to render the plot grid?

# Load data
WHITE = np.loadtxt("KERNEL-WHITE-05a.dat") # l, b_l
GAUSS = np.loadtxt("KERNEL-GAUSS-800.dat") # l, b_l
SSG21 = np.loadtxt("KERNEL-SSG21-800.dat") # l, b_l
SSG42 = np.loadtxt("KERNEL-SSG42-800.dat") # l, b_l
SSG84 = np.loadtxt("KERNEL-SSG84-800.dat") # l, b_l

# scale kernel to uK
WHITE[:,1] *= 1.0e-6

# Create the plots
for width in [18., 12., 8.8]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
    
    # beam functions
    plt.plot(WHITE[:,0], WHITE[:,1], "#B0C4DE", linewidth=3.0*width/8.8, label=r"\texttt{WHITE} kernel")
    plt.plot(SSG21[:,0], SSG21[:,1], "r", linewidth=1.5*width/8.8, label=r"\texttt{SSG21} kernel")
    plt.plot(SSG42[:,0], SSG42[:,1], "g", linewidth=1.5*width/8.8, label=r"\texttt{SSG42} kernel")
    plt.plot(SSG84[:,0], SSG84[:,1], "b", linewidth=1.5*width/8.8, label=r"\texttt{SSG84} kernel")
    plt.plot(GAUSS[:,0], GAUSS[:,1], "k", linewidth=1.5*width/8.8, label=r"\texttt{GAUSS} kernel")
    
    # x axis
    plt.hlines(0, 0, 4000, linewidth=0.5)
    
    # legend
    if (fill or width < 15):
        leg = plt.legend(frameon=True)
        # remove box around legend
        leg.get_frame().set_edgecolor('none')
        leg.get_frame().set_alpha(.8)
    else:
        #bbox = (0.0, 0.13, 1.0, 0.47) if width < 15.0 else (0.0, 0.0, 1.0, 1.0)
        #plt.legend(frameon=False, loc='upper right', bbox_to_anchor=bbox)
        #plt.legend(frameon=False, ncol=(2 if width < 15 else 1))
        plt.legend(frameon=False)
    
    # labels
    plt.xlabel(r"$\ell$"); plt.ylabel(r"$b_\ell$"); plt.title(r"FIR filter kernels (at $800^{\scriptstyle\prime}$ FWHM)")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
    
    # reduce ticks for small figures
    #if width < 10:
    #    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # grid
    if (grid):
        plt.grid(True, which="major", axis="both")
    
    # axes limits
    plt.ylim([-0.15, 1.05]); plt.xlim([0, 100]);
    
    # reduce white space around figure
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")
    
    # save to pdf with right bounding box
    plt.savefig("./peaks-beam-firs-800.%dmm.pdf" % int(width*10), bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
