module imageio ! common FITS image and matrix I/O operations

use healpix_types

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

public :: image2fits

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _dp; end interface

GENERIC(image2fits)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

! single precision
#define XP SP
#define VARIANT(name) name ## _sp
#include 'imageio.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _dp
#include 'imageio.fin'

end module
