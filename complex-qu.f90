module complex_qu ! complex representation of linearly polarized maps Z = Q + iU

! complex array slicing is implemented via pointers
use, intrinsic :: ISO_C_BINDING

! HEALPix modules
use healpix_types
use pix_tools
use udgrade_nr

implicit none

public :: convert_ring2nest, convert_nest2ring, udgrade_nest

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _zs, name ## _zd; end interface

GENERIC(convert_ring2nest)
GENERIC(convert_nest2ring)
GENERIC(udgrade_nest)

contains

! single precision
#define XP SP
#define VARIANT(name) name ## _zs
#include 'complex-qu.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _zd
#include 'complex-qu.fin'

end module
