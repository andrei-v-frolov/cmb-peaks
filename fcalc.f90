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
        case ('valid'); where (valid) Mout = 1.0
        case ('invalid'); where (.not. valid) Mout = 1.0
        case ('mask'); Mout = M1*M2; where (M2 == 0.0) Mout = 1.0/0.0
        case ('unmask'); Mout = M1/M2; where (M2 == 0.0) Mout = 1.0/0.0
        case ('inpaint'); call inpaint(M1, M2, Mout)
        
        ! unknown operator
        case default; call abort(trim(op) // ": operation not supported")
end select

! write output map
call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator='FCALC', version='$Revision$')
call output_map(Mout, header, '!'//fout)

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ...
subroutine inpaint(map, mask, mout)
        real(SP), dimension(0:n) :: map, mask, mout
        real(DP), dimension(0:n) :: M, L
        
        integer i, j, k, nn(9,0:n)
        real(DP) w(9,0:n)
        
        ! prepare stencil for the pixels that need inpainting
        k = 0; do i = 0,n
                M(i) = map(i)*mask(i); if (mask(i) /= 0.0) cycle
                call stencil(nside, ord, i, nn(:,k), bartjan=w(:,k))
                k = k+1
        end do
        
        ! iterate diffusion steps
        do j = 1,100
                forall (i=0:k-1) L(i) = sum(w(:,i)*M(nn(:,i)))
                M(nn(1,0:k-1)) = L(0:k-1)
        end do
        
        ! save the result
        mout = M
end subroutine inpaint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! return nearest neighbours list (in arbitrary map ordering)
subroutine neighbours(nside, order, i, nn, k)
        integer nside, order, i, j, k, nn(8)
        
        ! ordering convention literals
        integer, parameter :: RING = 1, NEST = 2
        
        ! get HEALPix neighbour list
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

! gnomonic coordinates (X,Y) of a pixel p around origin i (in arbitrary map ordering)
subroutine pix2gno(nside, order, i, p, XY)
        integer nside, order, i, p; real(DP) XY(2)
        real(DP) theta, phi, U(3), V(3), W(3), X(3), Y(3)
        
        ! ordering convention literals
        integer, parameter :: RING = 1, NEST = 2
        
        ! convert pixels to coordinates
        select case(order)
                case(RING)
                        call pix2vec_ring(nside, i, U)
                        call pix2vec_ring(nside, p, V)
                        call pix2ang_ring(nside, i, theta, phi)
                case(NEST)
                        call pix2vec_nest(nside, i, U)
                        call pix2vec_nest(nside, p, V)
                        call pix2ang_nest(nside, i, theta, phi)
                case default
                        call abort(": ordering not supported")
        end select
        
        ! project onto a tangent plane through origin
        W = V/sum(U*V) - U
        
        ! project onto local orthonormal basis in (theta,phi) directions
        X = (/ cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta) /)
        Y = (/ -sin(phi), cos(phi), 0.0 /)
        
        XY = (/ sum(W*X), sum(W*Y) /)
end subroutine pix2gno

! return differential operator stencil (to be applied to nearest neighbours)
subroutine stencil(nside, order, i, nn, count, bartjan)
        integer nside, order, i, j, k, nn(9)
        integer, intent(out), optional :: count
        
        ! supported stencils - see below for description
        real(DP), dimension(9), intent(out), optional :: bartjan
        
        ! nearest neighbour list
        nn(1) = i; call neighbours(nside, order, i, nn(2:), k)
        if (k < 8) nn(9) = i; if (present(count)) count = k+1
        
        ! nearest neighbour average used by Bartjan van Tent et. al.
        if (present(bartjan)) then; bartjan = 0.0; bartjan(2:k+1) = 1.0/k; end if
end subroutine stencil


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
