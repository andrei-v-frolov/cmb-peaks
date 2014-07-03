###############################################################################
# Gaussian random field peak statistics
###############################################################################

from math import *
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


###############################################################################
# statistics primitives
###############################################################################

def CDF(x, gamma, sigma=1.0, alpha=0.0):
    """Calculate fraction of Gaussian random field peaks (max+min) below threshold nu"""
    
    nu = x/sigma - alpha; mu = nu/sqrt(2.0-gamma**2/0.75)
    npk_above = sqrt(1.5/pi) * gamma**2 * nu*np.exp(-nu**2/2.0) + np.vectorize(erfc)(mu)/2.0
    npk_below = 1.0 - npk_above
    
    return npk_below

def lmoments(x, F):
    """Calculate lower L-moments of the distribution described by (tabulated) CDF F(x)""" 
    
    # weights are Legendre polynomials of percentile value
    u = 2.0*F - 1.0; v = 1.5*u**2 - 0.5; w = u*(2.5*u**2 - 1.5)
    
    return np.trapz(x, F), np.trapz(u*x, F), np.trapz(v*x, F), np.trapz(w*x, F)

def lkurtosis(gamma, n=1000, xmax=6.0):
    """Calculate L-kurtosis of a Gaussian peak distribution with parameter gamma"""
    
    x = np.linspace(-xmax, xmax, n); F = CDF(x, gamma)
    L1, L2, L3, L4 = lmoments(x, F)
    
    return gamma, L4/L2, L2-L4

# gamma and sigma estimators are based on L-kurtosis LUT
gamma_lut = np.vectorize(lkurtosis)(np.linspace(1.0, 0.0, 16))
gamma_est = interp1d(gamma_lut[1], gamma_lut[0], kind='cubic')
sigma_est = interp1d(gamma_lut[1], gamma_lut[2], kind='cubic')

def estimate(x, F):
    """Estimate parameters of a peak distribution F(x) from its L-moments"""
    
    L1, L2, L3, L4 = lmoments(x, F)
    
    gamma = gamma_est(L4/L2).item()
    sigma = (L2-L4)/sigma_est(L4/L2)
    alpha = (L1-L3)/sigma
    
    return gamma, sigma, alpha

def cdf_fit(x, F, n=32):
    """Determine parameters of a peak distribution F(x) by fitting Gaussian peak CDF"""
    
    # initial guess is calculated from L-moments
    fit = list(estimate(x, F)); s = len(x)/n + 1
    
    # refine the fit using progressively more data
    while (s > 0):
        fit, fitCovariances = curve_fit(CDF, x[::s], F[::s], fit); s >>= 1
    
    return fit, fitCovariances

def marginalize(func, p, cov, order=4):
    """Marginalize function of parameters with specified covariance matrix"""
    
    from itertools import product
    
    # integrate variance using Gauss-Hermite quadrature
    x, u = np.polynomial.hermite.hermgauss(order); v = np.log(u)
    
    # quadrature is evaluated on Cholesky-diagonalized variable z
    A = np.linalg.cholesky(cov)
    
    # best fit value and variance accumulator
    f = func(p); df = f - f
    
    # quadrature accumulation loop
    for i in product(range(order), repeat=len(p)):
        z = map(lambda k: x[k], i); dp = np.dot(A,z)
        q = map(lambda k: v[k], i);  w = exp(sum(q))
        
        df += w * (func(p+dp) - f)**2
    
    return f, np.sqrt(df)