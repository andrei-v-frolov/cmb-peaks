module almtools ! complex representation of linearly polarized maps Z = Q + iU

! complex array slicing is implemented via pointers
use, intrinsic :: ISO_C_BINDING

! HEALPix modules
use mapio
use alm_tools
use udgrade_nr

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HEALPix routine wrappers, covariant derivatives, in single and double precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _dp; end interface

GENERIC(pack_alms)
GENERIC(unpack_alms)
GENERIC(randomize_alms)
GENERIC(xrandomize_alms)
GENERIC(randomize)
GENERIC(alm2map_magnetic)
GENERIC(alm2map_covariant)
GENERIC(mask2spins_ring)

public :: pack_alms, unpack_alms
public :: randomize_alms, xrandomize_alms, randomize
public :: alm2map_magnetic, alm2map_covariant, mask2spins_ring

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HEALPix routine wrappers, complex QU, in single and double precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _zs, name ## _zd; end interface

GENERIC(convert_ring2nest)
GENERIC(convert_nest2ring)
GENERIC(udgrade_nest)

public :: convert_ring2nest, convert_nest2ring, udgrade_nest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! alm transform wrappers, real and complex QU, in single and double precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _dp, name ## _zs, name ## _zd; end interface

GENERIC(map2alm_pure)
GENERIC(rotate_qu2eb_pure)
GENERIC(rotate_qu2eb)
GENERIC(rotate_eb2qu)
GENERIC(project_qu)
GENERIC(randomize_qu)
GENERIC(xrandomize_qu)
GENERIC(krylov_qu)
GENERIC(purify_qu)

public :: map2alm_pure, rotate_qu2eb_pure, rotate_qu2eb, rotate_eb2qu, project_qu, randomize_qu, xrandomize_qu, krylov_qu, purify_qu

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! real QU maps
#define MAPTYPE(XP,DIM) real(XP), dimension(DIM,1:2)
#define LOAD(MAP,SRC) MAP = SRC
#define COPY(MAP,DST) DST = MAP

! single precision
#define XP SP
#define VARIANT(name) name ## _sp
#include 'almtools.fin'
#include 'maptools.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _dp
#include 'almtools.fin'
#include 'maptools.fin'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! complex QU maps
#define MAPTYPE(XP,DIM) complex(XP), dimension(DIM)
#define LOAD(MAP,SRC) MAP(:,1) = real(SRC); MAP(:,2) = imag(SRC)
#define COPY(MAP,DST) DST = cmplx(MAP(:,1), MAP(:,2))

! single precision
#define XP SP
#define VARIANT(name) name ## _zs
#include 'almtools.fin'
#include 'complex-qu.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _zd
#include 'almtools.fin'
#include 'complex-qu.fin'

end module
