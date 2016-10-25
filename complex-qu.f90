module complex_qu ! complex representation of linearly polarized maps Z = Q + iU

! complex array slicing is implemented via pointers
use, intrinsic :: ISO_C_BINDING

! HEALPix modules
use healpix_types
use pix_tools

implicit none

public :: convert_ring2nest

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _zs, name ## _zd; end interface

GENERIC(convert_ring2nest)

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
