module mapio ! common map I/O operations

use healpix_types
use fitstools
use head_fits
use pix_tools

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter, public :: IO = SP			! default I/O precision
integer, parameter, public :: RING = 1, NEST = 2	! ordering literals
integer, parameter, public :: CART = 1, HLPX = 2	! vector components
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! add vector field metadata
subroutine add_vector_card(header, vec)
	character(len=80) :: header(:); integer vec
	
	if (vec >= 0) call add_card(header, 'VECTOR', (vec > 0), 'Vector field included (True/False)')
	if (vec == CART) call add_card(header, 'VFRAME', 'CARTESIAN', 'Vector components frame')
	if (vec == HLPX) call add_card(header, 'VFRAME', 'HEALPIX',   'Vector components frame')
end subroutine

! get vector field metadata
function get_vector_card(header)
	character(len=80) :: header(:), frame
	integer get_vector_card, n; logical vector
	
	! no vector metadata
	get_vector_card = -1
	
	! is vector metadata present?
	call get_card(header, 'VECTOR', vector, count=n); if (n == 0) return
	if (.not. vector) then; get_vector_card = 0; return; end if
	
	! is coordinate frame specified?
	call get_card(header, 'VFRAME', frame,  count=n); if (n == 0) return
	
	select case (trim(frame))
		case ('CARTESIAN'); get_vector_card = CART
		case ('HEALPIX');   get_vector_card = HLPX
		case default; call abort(trim(frame) // ": coordinate frame not supported")
	end select
end function

end module
