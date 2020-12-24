#!/usr/bin/env python
# Plot pseudo-Cl spectra
# usage: plot-spectra <cls.fits> [output.pdf] [FWHM] [powerlaw.fits]

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

# FITS I/O libraries
try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits


###############################################################################
# plot style options
###############################################################################

width = 16.0	# figure width, in centimeters
aspect = 6/4.	# figure aspect ratio (default is 4:3)
fill = False	# produce transparent plots if false
grid = True	# do we want to render the plot grid?
legend = True	# do we want to include the legend?


###############################################################################
# parse arguments
###############################################################################
assert len(argv) > 1, 'usage: plot-spectra <cls.fits> [output.pdf] [FWHM] [powerlaw.fits]'

data = pyfits.open(argv[1])[1].data
file = argv[2] if len(argv) > 2 else None
fwhm = float(argv[3]) if len(argv) > 3 else 0.0
fits = argv[4] if len(argv) > 4 else None

l = np.arange(0,data.shape[0])

II = data['TEMPERATURE']
EE = data['GRADIENT']
BB = data['CURL']

IE = (data['TG'] + data['GT'])/2.0 if 'GT' in data.names else data['TG']
IB = (data['TC'] + data['CT'])/2.0 if 'CT' in data.names else data['TC']
EB = (data['GC'] + data['CG'])/2.0 if 'CG' in data.names else data['GC']

###############################################################################
# logarithmic binning
###############################################################################

def binned(l,Cl):
	"""Compute average ell, power, and variance in a bin"""
	P = np.sum((2*l+1)*Cl); N = np.sum(2*l+1); L = np.sum(l)
	return L/len(l), P/N, sqrt(np.sum((2*l+1)*(Cl-P/N)**2)/(N-1))

def logbin(l,Cl,bands=100):
	"""Bin spectrum data in logarithmically spaced bins"""
	n = len(l); delta = np.exp(np.log(l[-1]/max(l[0],2.0))/bands) - 1.0
	
	i = 0; data = []; bins = []; sigma = []
	
	while (i < n):
		s = min(2*int(l[i]*delta) + 1, n-i)
		b,v,w = binned(l[i:i+s],Cl[i:i+s])
		bins.append(b); data.append(v); sigma.append(w)
		i += s
	
	return np.array(bins),np.array(data),np.array(sigma)

###############################################################################
# fitting data
###############################################################################

from warnings import warn
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

# effective beam (20' FWHM Gaussian/Pixel window at n=2048)
BFWHM = fwhm * pi/(180*60) # 20 arcmin
P2048 = 1.21333461333613 * pi/(180*60)
sigma = (BFWHM**2-P2048**2)/(8.0*log(2.0))

def powerlaw(x, a, alpha, sigma=sigma):
	return a * np.power(x, alpha-1)/(x+1) * np.exp(-sigma*x*(x+1))

def bestfit(l, XX, lmin=20, lmax=200):
	# approximate first guess
	x,y,s = logbin(l[lmin:lmax],XX[lmin:lmax],bands=40); mask = y > 0.0
	alpha,b,r,p,e = linregress(np.log(x[mask]), np.log(y[mask]))
	
	alpha += 2; a = exp(b)
	
	#def residual(p):
	#	return (y - np.vectorize(powerlaw)(x, p[0], p[1]))/1.0e-11
	#fit = least_squares(residual, np.array([a, alpha]), loss='arctan')
	#return fit.x[1], np.vectorize(powerlaw)(l, fit.x[0], fit.x[1])
	
	# curve fit might fail, fallback to approximate guess
	try:
		fit,cov = curve_fit(powerlaw, x, y, [a, alpha], sigma=s)
	except RuntimeError:
		warn("curve fit failed to converge, falling back to approximate values")
		fit = [a, alpha]
	
	return fit, np.vectorize(powerlaw)(l, fit[0], fit[1])

###############################################################################
# create the plots
###############################################################################

os.environ['matplotlib.backend'] = 'macosx' if file is None else 'pdf'

# Configure Matplotlib options
from pltconfig import *

fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/aspect)), frameon=fill)

k,v,w = binned(l[100:500],EB[100:500])
ymax = 20.0 * 10.0**ceil( log10(0.5 * max((II * l*(l+1)/(2.0*np.pi))[2:100])))
ymin = 0.05 * 10.0**floor(log10(2.0 * sqrt(v*v + w*w/10.0) * k*(k+1)/(2.0*np.pi)))

###############################################################################

current = plt.subplot(1,2,1); a = 0; b = -20

plt.loglog(l, II * l*(l+1)/(2.0*np.pi), 'r-', label='', alpha=0.2)
plt.loglog(l, EE * l*(l+1)/(2.0*np.pi), 'g-', label='', alpha=0.2)
plt.loglog(l, BB * l*(l+1)/(2.0*np.pi), 'b-', label='', alpha=0.2)

k,XX,_ = logbin(l[a:b],II[a:b]); plt.loglog(k, XX * k*(k+1)/(2.0*np.pi), 'r-', label='ii' if 'log' in argv[1] else 'II', alpha=0.5)
k,XX,_ = logbin(l[a:b],EE[a:b]); plt.loglog(k, XX * k*(k+1)/(2.0*np.pi), 'g-', label='ee' if 'log' in argv[1] else 'EE', alpha=0.5)
k,XX,_ = logbin(l[a:b],BB[a:b]); plt.loglog(k, XX * k*(k+1)/(2.0*np.pi), 'b-', label='bb' if 'log' in argv[1] else 'BB', alpha=0.5)

IIfit,XX = bestfit(l, II); plt.loglog(l, XX * l*(l+1)/(2.0*np.pi), 'r--', label=(r'$\alpha = %3.2f$' % (IIfit[1]-2.0)))
EEfit,XX = bestfit(l, EE); plt.loglog(l, XX * l*(l+1)/(2.0*np.pi), 'g--', label=(r'$\alpha = %3.2f$' % (EEfit[1]-2.0)))
BBfit,XX = bestfit(l, BB); plt.loglog(l, XX * l*(l+1)/(2.0*np.pi), 'b--', label=(r'$\alpha = %3.2f$' % (BBfit[1]-2.0)))

if grid: plt.grid(True, which="both", axis="both", linestyle='-', linewidth=0.2)
if legend: plt.legend(loc='upper right', ncol=2, frameon=fill or grid)

plt.ylabel(r'$\ell(\ell+1)\ C_\ell^{XX}/2\pi$')
plt.xlabel(r'$\ell$')

plt.xticks([2,10,100,1000], ['2','10','100','1000'])
plt.xlim(2,1000); plt.ylim(ymin,ymax)

###############################################################################

current = plt.subplot(1,2,2)
current.yaxis.tick_right()
current.yaxis.set_label_position("right")

plt.loglog(l, abs(IE) * l*(l+1)/(2.0*np.pi), '-', label='', color="Orange", alpha=0.2)
plt.loglog(l, abs(IB) * l*(l+1)/(2.0*np.pi), '-', label='', color="DarkCyan", alpha=0.2)
plt.loglog(l, abs(EB) * l*(l+1)/(2.0*np.pi), '-', label='', color="Purple", alpha=0.2)

k,XX,_ = logbin(l[a:b],IE[a:b]); plt.loglog(k, abs(XX) * k*(k+1)/(2.0*np.pi), '-', label='ie' if 'log' in argv[1] else 'IE', color="Orange", alpha=0.5)
k,XX,_ = logbin(l[a:b],IB[a:b]); plt.loglog(k, abs(XX) * k*(k+1)/(2.0*np.pi), '-', label='ib' if 'log' in argv[1] else 'IB', color="DarkCyan", alpha=0.5)
k,XX,_ = logbin(l[a:b],EB[a:b]); plt.loglog(k, abs(XX) * k*(k+1)/(2.0*np.pi), '-', label='eb' if 'log' in argv[1] else 'EB', color="Purple", alpha=0.5)

IEfit,XY = bestfit(l, abs(IE)); plt.loglog(l, XY * l*(l+1)/(2.0*np.pi), '--', label=(r'$\alpha = %3.2f$' % (IEfit[1]-2.0)), color="Orange")
IBfit,XY = bestfit(l, abs(IB)); plt.loglog(l, XY * l*(l+1)/(2.0*np.pi), '--', label=(r'$\alpha = %3.2f$' % (IBfit[1]-2.0)), color="DarkCyan")
EBfit,XY = bestfit(l, abs(EB)); plt.loglog(l, XY * l*(l+1)/(2.0*np.pi), '--', label=(r'$\alpha = %3.2f$' % (EBfit[1]-2.0)), color="Purple")

if grid: plt.grid(True, which="both", axis="both", linestyle='-', linewidth=0.2)
if legend: plt.legend(loc='upper right', ncol=2, frameon=fill or grid)

plt.ylabel(r'$\ell(\ell+1)\ |C_\ell^{XY}|/2\pi$')
plt.xlabel(r'$\ell$')

plt.xticks([2,10,100,1000], ['2','10','100','1000'])
plt.xlim(2,1000); plt.ylim(ymin,ymax)

###############################################################################

plt.tight_layout()
plt.show()

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=not(fill))


###############################################################################
# output fitted spectra
###############################################################################

if fits is None: sys.exit()

def column(name, fit, lmin=0, lmax=3000, format='E15.7', units='unknown^2', replace=[]):
	data = np.vectorize(powerlaw)(np.arange(0,lmax+1), fit[0], fit[1], sigma=0.0)
	if len(replace) > 0: data[0:len(replace)] = replace
	if lmin > 0: data[0:lmin] = 0.0
	return pyfits.Column(name, format, units, array=data)

columns = [
	column('TEMPERATURE', IIfit, replace=II[0:1]),
	column('GRADIENT', EEfit, 2),
	column('CURL', BBfit, 2),
	column('TG', IEfit, 2),
	column('TC', IBfit, 2),
	column('GC', EBfit, 2)
]

hdu = pyfits.PrimaryHDU(pyfits.Header())
tbl = pyfits.TableHDU.from_columns(columns)
pyfits.HDUList([hdu,tbl]).writeto(fits)
