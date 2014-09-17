###############################################################################
# import Planck-specific colormaps
###############################################################################

# configure import path
import os; here = os.path.split(os.path.realpath(__file__))[0]

# import libraries
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.cm import register_cmap

# colormap loader
def load_cmap(name, file, bad='gray', over='white', under='white'):
    """Load and register colormap from an RGB data file"""
    
    cmap = ListedColormap(np.loadtxt(file)/255.0, name=name)
    cmap.set_bad(bad); cmap.set_over(over); cmap.set_under(under)
    
    register_cmap(name, cmap)
    
    return cmap

# default HEALPix and Planck colormaps
healpix_cmap = load_cmap('healpix', os.path.join(here, 'colormaps/HEALPix_CMB.rgb'))
planck_cmap  = load_cmap('parchment', os.path.join(here, 'colormaps/Planck_Parchment.rgb'))
