#!/usr/bin/env python
# Plot pre-whitener kernels; usage: beam-whiten [width list (in cm)]

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

widths = map(float, sys.argv[1].split(',')) if (len(sys.argv) > 1) else [18.0, 12.0, 8.8]

# Load data
WHITE = np.loadtxt("KERNEL-WHITE-05a.dat") # l, b_l

# scale kernel to uK
WHITE[:,1] *= 1.0e-6

# Create the plots
for width in widths:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
    
    # beam functions
    plt.plot(WHITE[:,0], WHITE[:,1], "#B0C4DE", linewidth=3.0*width/8.8, label=r"\texttt{WHITE} kernel")
    
    # x axis
    plt.hlines(0, 0, 4000, linewidth=0.5)
    
    # legend
    if (fill):
        leg = plt.legend(frameon=True)
        # remove box around legend
        leg.get_frame().set_edgecolor('none')
        leg.get_frame().set_alpha(.8)
    else:
        plt.legend(frameon=False)
    
    # labels
    plt.xlabel(r"$\ell$"); plt.ylabel(r"$b_\ell$"); plt.title(r"Pre-whiten filter kernel, full $\ell$ range")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
    
    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # grid
    if (grid):
        plt.grid(True, which="major", axis="both")
    
    # axes limits
    plt.ylim([0.0, 40.0]); plt.xlim([0, 4000]);
    
    # reduce white space around figure
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")
    
    # save to pdf with right bounding box
    base = 'peaks-beam-whiten'
    pdf = "%s.%dmm.pdf" % (base,int(width*10)) if len(widths) > 1 else base+".pdf"
    plt.savefig(pdf, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
