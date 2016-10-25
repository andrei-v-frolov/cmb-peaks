module mapio ! common map I/O operations

use healpix_types
use fitstools
use head_fits
use pix_tools

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter, public :: IO = SP			! default I/O precision
integer, parameter, public :: RING = 1, NEST = 2	! ordering literals
logical, parameter, public :: verbose = .true.		! diagnostic output

public :: read_map, write_map

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _dp; end interface

GENERIC(read_map)
GENERIC(write_map)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

! single precision
#define XP SP
#define VARIANT(name) name ## _sp
#include 'mapio.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _dp
#include 'mapio.fin'

! output warning message to stderr
subroutine warning(msg)
	character(*) msg; logical, parameter :: verbose = .true.
	
	if (verbose) write (0,'(a)') "warning: " // msg
end subroutine warning

end module
