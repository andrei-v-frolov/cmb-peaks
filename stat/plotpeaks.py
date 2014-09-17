#!/usr/bin/env python
# Plot peak distribution and deviation from Gaussian peak statistics
# usage: peakplot KERNEL FWHM [input datafile] [plot basename] [width list (in cm)]

###############################################################################
# import libraries
###############################################################################

# configure import path
import os, sys
here = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(here, '../libs'))

# import libraries
from sys import argv, stdin
from peakstats import *
from spherical import *

# plot style options
aspect = 4/3.   # figure aspect ratio (default is 4:3)
fill = False    # produce transparent plots if false
grid = False    # do we want to render the plot grid?
ptlimit = 513   # point limit for scatter plots (decimate if necessary)

# dipole modulation direction (set to None for full sky only)
dipole = iau2vec(227,-15); # dipole = None

# parse arguments
kernel = argv[1]
fwhm = float(argv[2])
file = argv[3] if (len(argv) > 3) else "L1-%s-%03i" % (kernel, fwhm/2) + ".dat"
base = argv[4] if (len(argv) > 4) else "L1-%s-%03i" % (kernel, fwhm/2)
widths = map(float, argv[5].split(',')) if (len(argv) > 5) else [18.0, 12.0, 8.8]


###############################################################################
# import peak data
###############################################################################

# peak data format: theta, phi, value, kind (-1 => min, +1 => max);
peaks = np.loadtxt(file if file != '-' else stdin)

# form peak CDF
x, f, n = makecdf(peaks[:,2])

# fit Gaussian random peak distribution
fit, cov = cdf_fit(x,f)
gamma, sigma, alpha = fit


###############################################################################
# North/South asymmetry
###############################################################################

north = []
south = []

if (dipole != None):
    for i in range(n):
        q = np.dot(ang2vec(peaks[i,0], peaks[i,1]), dipole)
        
        if (q < -1.0/3.0):
            north.append(peaks[i,2])
        if (q > +1.0/3.0):
            south.append(peaks[i,2])


###############################################################################
# format plot labels
###############################################################################

def signlbl(loglikelihood):
    s, ds = loglikelihood
    return r"\large\textbf{significance 1:%.1f$^{%+.1f}_{%+.1f}$}" % (exp(-s), exp(-s+ds)-exp(-s), exp(-s-ds)-exp(-s))

fitlbl =r"Gaussian peak fit, $\gamma = %.2f$" % gamma
datlbl = r"\texttt{%s} filter at $%.0f^{\scriptstyle\prime}$ {\small FWHM}" % (kernel, fwhm)
count = r"\underline{\Large\bf %i peaks}" % n

coldspot = r"coldest at $(%+.2f%+.2f)\sigma$" % (alpha, x[0]/sigma-alpha)
coldsign = marginalize(lambda p: log(1.0 - (1.0-CDF(x[0], p[0], p[1], p[2]))**n), fit, cov)
coldspot += "\n" + signlbl(coldsign)

hotspot = r"hottest at $(%+.2f%+.2f)\sigma$" % (alpha, x[-1]/sigma-alpha)
hotsign = marginalize(lambda p: log(1.0 - CDF(x[-1], p[0], p[1], p[2])**n), fit, cov)
hotspot += "\n" + signlbl(hotsign)


###############################################################################
# create the plots
###############################################################################

# Configure Matplotlib options
from pltconfig import *
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

def plot_setup(ax, xlabel='', ylabel='', xlim=[-6,6], ylim=[0,1], legend=''):
    """Setup common plot style parameters"""
    
    # legend
    if (legend):
        if (fill):
            leg = plt.legend(loc=legend, frameon=True)
            # remove box around legend
            leg.get_frame().set_edgecolor('none')
            leg.get_frame().set_alpha(.8)
        else:
            plt.legend(loc=legend, frameon=False)
    
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

def cdfpts(data, color='r', marker='+', size=1.0, label='', zorder=5):
    """Plot data points constituting CDF, decimating if necessary"""
    
    # form peak CDF from value list
    x, f, n = makecdf(data)
    
    # decimate distribution
    if (ptlimit and n > 3*ptlimit/2):
        x,f = decimate(x,f,n=ptlimit)
    
    # plot data
    plt.scatter(x/sigma, f, 24*width/8.8*size, color=color, marker=marker, linewidth=0.5, label=label, zorder=zorder)

def kstest(data, color='r', marker='+', size=1.0, label='', zorder=5):
    """Plot Kolmogorov-Smirnov deviation for a specified data set"""
    
    # form peak CDF from value list
    x, f, n = makecdf(data)
    
    # decimate distribution
    if (ptlimit and n > 3*ptlimit/2):
        x,f = decimate(x,f,n=ptlimit)
    
    # Kolmogorov-Smirnov deviation
    K = sqrt(n)+0.12+0.11/sqrt(n)
    
    # evaluate best fit CDF and fit variance
    y, dy = marginalize(lambda p: CDF(x, p[0], p[1], p[2]), fit, cov)
    
    # plot deviation from Gaussian peak CDF
    for i in range(3):
        plt.fill_between(x/sigma, K*(f-y-(i+1)*dy), K*(f-y+(i+1)*dy), edgecolor='none', color=color, alpha=0.1, zorder=zorder-(i+1))
    
    if (marker):
        plt.scatter(x/sigma, K*(f-y), 24*width/8.8*size, color=color, marker=marker, linewidth=0.5, label=label, zorder=zorder)
    #plt.errorbar(x/sigma, K*(f-y), yerr=3*K*dy, color="r", fmt="o", linewidth=0.5, label="Data")

def plot_ks_panel(plot, xlim=[-6,6]):
    """Assemble Kolmogorov-Smirnov deviation panel plots"""
    
    plt.title("Kolmogorov deviation from Gaussian peak CDF")
    plt.hlines(0, xlim[0], xlim[1], linewidth=0.5, zorder=0) # x axis
    
    # Kolmogorov-Smirnov confidence levels
    for i in [0.3989178859, 0.5137003373, 0.7170542898]:
        plt.fill_between(xlim, -i, i, edgecolor='none', color="darkgrey", alpha=0.2, zorder=-5)
    
    if (dipole != None):
        kstest(north, color='b', marker='x', size=0.2, label="northern cap")
        kstest(south, color='g', marker='x', size=0.2, label="southern cap")
        kstest(x, marker='+', size=0.5, label="entire sky")
        
        plot_setup(plot, ylim=[-1.3,1.3])
        plt.legend(loc='upper center', ncol=3, frameon=False)
    else:
        kstest(x, marker='+', size=0.5)
        
        plot_setup(plot, ylim=[-1,1])
        plt.yticks([-0.5,0.0,0.5,1.0])

def plot_cdf_panel(plot, xlim=[-6,6], n=1024):
    """Assemble peak CDF panel plots"""
    
    from matplotlib.ticker import FormatStrFormatter
    
    # plot best-fitting Gaussian CDF
    nu = np.linspace(xlim[0], xlim[1], n)
    y, dy = marginalize(lambda p: CDF(sigma*nu, p[0], p[1], p[2]), fit, cov)
    
    for i in range(3):
        plt.fill_between(nu, y-(i+1)*dy, y+(i+1)*dy, edgecolor='none', color="g", alpha=0.1, zorder=-(i+1))
    plt.plot(nu, y, "g", linewidth=1.5*width/8.8, zorder=0, label=fitlbl)
    
    # plot actual peak distribution
    cdfpts(x, label=datlbl)
    
    # annotate with pretty labels
    if (width > 15):
        plt.annotate(coldspot, [xlim[0]+0.1,0.05], va='bottom', ha='left', linespacing=1.5)
        plt.annotate(hotspot,  [xlim[1]-0.1,0.95], va='top', ha='right', linespacing=1.5)
        plt.annotate(count, [xlim[1]-0.3,0.05], va='bottom', ha='right', linespacing=1.5)
    else:
        plt.annotate((r"\small " if width < 10 else '') + "coldest", xy = (x[0]/sigma, 0.01), xytext = (0, 20),
            textcoords = 'offset points', ha = 'center', va = 'bottom', color='white',
            bbox = dict(boxstyle = 'round,pad=0.3', ec='none', fc = 'blue', alpha = 0.2),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        plt.annotate((r"\small " if width < 10 else '') + "hottest", xy = (x[-1]/sigma, 0.99), xytext = (0, -20),
            textcoords = 'offset points', ha = 'center', va = 'top', color='white',
            bbox = dict(boxstyle = 'round,pad=0.3', ec='none', fc = 'red', alpha = 0.2),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        plt.annotate(count +'\n' + (hotspot+"\n" if width > 10 else '')+coldspot, [xlim[1]-0.3,0.05], va='bottom', ha='right', linespacing=1.5)
    
    plot_setup(plot, xlim=xlim, legend='upper left')
    plot.xaxis.set_major_formatter(FormatStrFormatter('$%d\sigma$'))

# Create the plots
for width in widths:
    layout = gridspec.GridSpec(2, 1, height_ratios=[5, 8])
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)
    
    # create plots
    cdf = plt.subplot(layout[1]); plot_cdf_panel(cdf)
    ks  = plt.subplot(layout[0], sharex=cdf); plot_ks_panel(ks)
    
    # x axes is shared, so we don't need ticks except in the last plot
    plt.setp([a.get_xticklabels() for a in fig.axes[1:]], visible=False)
    
    # reduce white space around figure
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0)
    
    # save to pdf with right bounding box
    pdf = "%s.%dmm.pdf" % (base,int(width*10)) if len(widths) > 1 else base+".pdf"
    plt.savefig(pdf, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))
