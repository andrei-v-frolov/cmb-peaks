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

! maptools.fin
GENERIC(pack_alms)
GENERIC(unpack_alms)
GENERIC(lconvolution)
GENERIC(map2grad)
GENERIC(alm2map_covariant)
GENERIC(mask2spins_ring)

! magnetic.fin
GENERIC(magnetic2pqu)
GENERIC(pqu2magnetic)
GENERIC(magnetic_fit)
GENERIC(magnetic_mcmc)
GENERIC(alm2map_magnetic)

! wavelets.fin
GENERIC(einbeam)
GENERIC(sample_beam)
GENERIC(linterp1d)
GENERIC(interpolate)
GENERIC(huang)

! randomize.fin
GENERIC(randomize_alms)
GENERIC(xrandomize_alms)
GENERIC(randomize)

public :: pack_alms, unpack_alms, lconvolution, map2grad, alm2map_covariant, mask2spins_ring
public :: magnetic2pqu, pqu2magnetic, magnetic_fit, magnetic_mcmc, alm2map_magnetic
public :: einbeam, sample_beam, linterp1d, interpolate, huang
public :: randomize_alms, xrandomize_alms, randomize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HEALPix routine wrappers, complex QU, in single and double precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _zs, name ## _zd; end interface

! complex-qu.fin
GENERIC(convert_ring2nest)
GENERIC(convert_nest2ring)
GENERIC(udgrade_nest)

public :: convert_ring2nest, convert_nest2ring, udgrade_nest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! alm transform wrappers, real and complex QU, in single and double precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _dp, name ## _zs, name ## _zd; end interface

! almtools.fin
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
#include 'magnetic.fin'
#include 'wavelets.fin'
#include 'randomize.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _dp
#include 'almtools.fin'
#include 'maptools.fin'
#include 'magnetic.fin'
#include 'wavelets.fin'
#include 'randomize.fin'

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
