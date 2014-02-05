! $Id$
! Calculate a difference between two maps
! invoke: fdiff A.fits B.fits output.fits 

program fdiff

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: nmaps = 0, nside = 0, ord = 0
character(len=80) :: header(64), fin1, fin2, fout
real(SP), dimension(:,:), allocatable :: M1, M2, Mout


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
call getArgument(1, fin1)
call getArgument(2, fin2)
call getArgument(3, fout)

! read input maps
call read_map(fin1, M1, nside, nmaps, ord)
call read_map(fin2, M2, nside, nmaps, ord)

allocate(Mout, mold=M1)

Mout = M1/M2

! write output map
call write_minimal_header(header, 'MAP', nside=nside, order=ord)
call output_map(Mout, header, '!'//fout)

contains

! read map from FITS file, allocating storage if necessary
subroutine read_map(fin, M, nside, nmaps, ord)
        character(*) fin
        real(SP), allocatable :: M(:,:)
        integer nside, npix, nmaps, ord
        
        ! read header info
        character(len=80) :: header(64)
        integer hside, htot, hmaps, hord
        integer, parameter :: RING = 1, NEST = 2
        
        htot = getsize_fits(fin, nside=hside, nmaps=hmaps, ordering=hord)
        if (htot == -1) call abort(trim(fin) // ": file not found")
        
        ! check if map format agrees with requested one
        if (nside == 0) nside = hside; if (hside /= nside) call abort(trim(fin) // ": map resolution does not conform")
        if (nmaps == 0) nmaps = hmaps; if (hmaps  < nmaps) call abort(trim(fin) // ": too few channels in an input map")
        if (  ord == 0)   ord = hord;  if ( hord /= ord) call warning(trim(fin) // ": map ordering is being converted")
        
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

end
