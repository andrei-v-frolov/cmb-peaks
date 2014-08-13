! $Id$
! Calculate a pointwise binary operation on two maps
! invoke: fcalc A.fits 'x' B.fits ['=>'] output.fits
! x is a binary operator, see source code for complete list

program fcalc

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: nmaps = 0, nside = 0, ord = 0, n = 0
character(len=80) :: header(64), fin1, op, fin2, fout
real(SP), dimension(:,:), allocatable :: M1, M2, Mout
logical, dimension(:,:), allocatable :: valid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
call getArgument(1, fin1)
call getArgument(2, op  )
call getArgument(3, fin2)
call getArgument(4, fout)

! syntax sugar
if (fout .eq. '=>') call getArgument(5, fout)

! read input maps
call read_map(fin1, M1, nside, nmaps, ord)
call read_map(fin2, M2, nside, nmaps, ord)
n = nside2npix(nside)-1

! output storage
allocate(Mout, mold=M1); Mout = 0.0

! valid data mask
allocate(valid(0:n,nmaps))
valid = .not. (isnan(M1) .or. isnan(M2))

! apply operator
select case (op)
        ! arithmetic operators
        case ('+');  Mout = M1 + M2
        case ('-');  Mout = M1 - M2
        case ('*');  Mout = M1 * M2
        case ('/');  Mout = M1 / M2
        case ('//'); Mout = floor(M1/M2)
        case ('**'); Mout = M1 ** M2
        
        ! comparison operators
        case ('<');  where (M1 < M2) Mout = 1.0
        case ('>');  where (M1 > M2) Mout = 1.0
        case ('<='); where (M1 <= M2) Mout = 1.0
        case ('>='); where (M1 >= M2) Mout = 1.0
        case ('=','=='); where (M1 == M2) Mout = 1.0
        case ('!=','/=','<>'); where (M1 /= M2) Mout = 1.0
        
        ! projection operators
        case ('project on'); Mout = sum(M1*M2,valid)/sum(M2*M2,valid) * M2
        case ('orthogonal'); Mout = M1 - sum(M1*M2,valid)/sum(M2*M2,valid) * M2
        
        ! masking operators
        case ('mask'); Mout = M1*M2; where (M2 == 0.0) Mout = 1.0/0.0
        case ('valid'); where (valid) Mout = 1.0
        case ('inpaint'); call inpaint(M1, M2, Mout)
        
        ! unknown operator
        case default; call abort(trim(op) // ": operation not supported")
end select

! write output map
call write_minimal_header(header, 'MAP', nside=nside, order=ord)
call output_map(Mout, header, '!'//fout)

contains

! ...
subroutine inpaint(map, mask, mout)
        real(SP), dimension(0:n) :: map, mask, mout
        real(DP), dimension(0:n) :: M, L
        real(DP) stencil(9,0:n)
        integer i, k, nn(9,0:n)
        
        M = map*mask; stencil = 0.0
        
        do i = 0,n
                call neighbours(nside, i, nn(2:9,i), k, ord)
                nn(1,i) = i; if (mask(i) /= 0.0) k = 0
                
                stencil(1:k+1,i) = 1.0/(k+1)
        end do
        
        do k = 1,100
                forall (i=0:n) L(i) = sum(stencil(:,i)*M(nn(:,i)))
                where (mask == 0.0) M = L
        end do
        
        mout = M
end subroutine inpaint

! return nearest neighbours list in arbitrary ordering
subroutine neighbours(nside, i, nn, k, order)
        integer, parameter :: RING = 1, NEST = 2
        integer i, j, k, nn(8), nside, order
        
        select case(order)
                case(RING)
                        call ring2nest(nside, i, j)
                        call neighbours_nest(nside, j, nn, k)
                        do j = 1,k
                                call nest2ring(nside, nn(j), nn(j))
                        end do
                case(NEST)
                        call neighbours_nest(nside, i, nn, k)
                case default
                        call abort(": ordering not supported")
        end select
end subroutine neighbours

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
