module mapio ! common map I/O operations

use healpix_types
use fitstools
use head_fits
use pix_tools

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter, public :: IO = SP                   ! default I/O precision
integer, parameter, public :: RING = 1, NEST = 2        ! ordering literals
logical, parameter, public :: verbose = .true.          ! diagnostic output

public :: read_map, write_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

! write map to FITS file, creating minimal header
subroutine write_map(fout, M, nside, ord, creator, version)
        character(*) fout, creator, version
        real(IO) M(:,:); integer nside, ord
        optional creator, version

        ! FITS header
        character(len=80) :: header(64)
        
        ! mode of operation
	integer mode; mode = 0
	if (present(creator)) mode = mode + 2
	if (present(version)) mode = mode + 1
	
	! write output header
	select case(mode)
		case(3); call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator=creator, version=version)
		case(2); call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator=creator)
		case(1); call write_minimal_header(header, 'MAP', nside=nside, order=ord, version=version)
		case default; call write_minimal_header(header, 'MAP', nside=nside, order=ord)
	end select
	
	! write output map
	call output_map(M, header, '!'//fout)
end subroutine write_map

! read map from FITS file, allocating storage if necessary
subroutine read_map(fin, M, nside, nmaps, ord)
        character(*) fin
        real(IO), allocatable :: M(:,:)
        integer nside, npix, nmaps, ord
        
        ! read header info
        character(len=80) :: header(64)
        integer hside, htot, hmaps, hord
        
        htot = getsize_fits(fin, nside=hside, nmaps=hmaps, ordering=hord)
        if (htot == -1) call abort(trim(fin) // ": file not found")
        
        ! check if map format agrees with requested one
        if (nside == 0) nside = hside; if (hside /= nside) call abort(trim(fin) // ": map resolution does not conform")
        if (nmaps == 0) nmaps = hmaps; if (hmaps  < nmaps) call abort(trim(fin) // ": too few channels in an input map")
                                       if (hmaps  > nmaps) call warning(trim(fin) // ": ignoring extra channels")
        if (  ord == 0)   ord = hord;  if ( hord /= ord)   call warning(trim(fin) // ": map ordering is being converted")
        
        ! allocate storage if needed
        npix = nside2npix(nside); if (.not. allocated(M)) allocate(M(npix,nmaps))
        if (size(M,1) /= npix .or. size(M,2) < nmaps) call abort(trim(fin) // ": unexpected storage array shape")
        
        ! read map data, converting order if needed
        call input_map(fin, M, npix, nmaps)
        if (hord == RING .and. ord == NEST) call convert_ring2nest(nside, M)
        if (hord == NEST .and. ord == RING) call convert_nest2ring(nside, M)
end subroutine read_map

! output warning message to stderr
subroutine warning(msg)
        character(*) msg; logical, parameter :: verbose = .true.
        
        if (verbose) write (0,'(a)') "warning: " // msg
end subroutine warning

end module