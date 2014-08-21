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

integer, parameter :: IO = SP               ! default I/O precision
integer, parameter :: RING = 1, NEST = 2    ! ordering literals

integer :: nmaps = 0, nside = 0, ord = 0, n = 0
character(len=80) :: header(64), fin1, op, fin2, fout
real(IO), dimension(:,:), allocatable :: M1, M2, Mout
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
        case ('inpaint'); call inpaint_mg(M1, M2, Mout)
        
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
        real(IO), dimension(0:n) :: map, mask, mout
        real(DP), dimension(0:n) :: U, V
        
        integer i, k, m, nn(9,0:n)
        real(DP) w(9,0:n)
        
        ! prepare stencil for the pixels that need inpainting
        k = 0; do i = 0,n
                U(i) = map(i)*mask(i); if (mask(i) /= 0.0) cycle
                call stencil(nside, ord, i, nn(:,k), bartjan=w(:,k))
                k = k+1
        end do; m = k-1
        
        ! iterate diffusion steps
        do i = 1,100
                forall (k=0:m) V(k) = sum(w(:,k)*U(nn(:,k)))
                forall (k=0:m) U(nn(1,k)) = V(k)
        end do
        
        ! save the result
        mout = U
end subroutine inpaint

! ...
subroutine inpaint_mg(map, mask, mout)
        real(IO), dimension(0:n) :: map, mask, mout
        real(DP), dimension(0:n) :: U, K, R
        
        integer i
        
        U = map*mask; K = mask
        
        select case(ord)
                case(RING)
                        call convert_ring2nest(nside, U)
                        call convert_ring2nest(nside, K)
                        
                        do i = 1,32
                            R = 0.0; call wstroke(nside, U, K, R, R)
                            write (*,*) "Total residual", sqrt(sum(R*R))
                        end do
                        
                        call convert_nest2ring(nside, U)
                        call convert_nest2ring(nside, R)
                case(NEST)
                        do i = 1,32
                            R = 0.0; call wstroke(nside, U, K, R, R)
                            write (*,*) "Total residual", sqrt(sum(R*R))
                        end do
                case default
                        call abort(": ordering not supported")
        end select
        
        mout = U
end subroutine inpaint_mg

! multigrid W-stroke: inpaints with L[map] = rhs where mask is not unity
! all maps are assumed to be in nested ordering for performance reasons
recursive subroutine wstroke(nside, map, mask, rhs, residual); use udgrade_nr
        integer i, k, m, n, nside, niter; real(DP) h2
        
        ! fine and coarse maps
        real(DP), dimension(0:12*nside**2-1) :: map, mask, rhs, tmp
        real(DP), dimension(0:12*nside**2-1), optional :: residual
        real(DP), dimension(0:3*nside**2-1) :: cmap, cmask, crhs
        intent(INOUT) map; intent(in) mask, rhs; intent(out) residual
        
        ! stencil operators
        integer,  dimension(9,0:12*nside**2-1) :: nn    ! nearest neighbours
        real(DP), dimension(9,0:12*nside**2-1) :: L     ! Laplacian stencil
        
        ! init pixel ranges
        n = nside2npix(nside)-1; m = 0; k = 0
        
        ! pixel grid spacing
        h2 = (pi/3.0)/nside**2
        
        ! Jacobi smoothing schedule
        niter = 64; if (nside > 8) niter = 16
        
        !write (*,*) "Entering w-stroke at nside=", nside
        
        ! stencils for the pixels that need inpainting
        do i = 0,n; if (mask(i) == 1.0) cycle
                call stencil(nside, NEST, i, nn(:,k), laplace=L(:,k)); k = k+1
        end do; m = k-1
        
        ! pre-smooth using Jacobi iteration
        do i = 1,niter
                forall (k=0:m) tmp(k) = (h2*rhs(nn(1,k)) - sum(L(2:9,k)*map(nn(2:9,k))))/L(1,k)
                forall (k=0:m) map(nn(1,k)) = tmp(k)
        end do
        
        ! solve coarse problem recursively (if we are not on the coarsest grid, that is)
        if (nside > 8) then
                ! downgrade residual
                cmap = 0.0; tmp = 0.0
                
                forall (k=0:m) tmp(nn(1,k)) = rhs(nn(1,k)) - sum(L(:,k)*map(nn(:,k)))/h2
                
                call udgrade_nest(mask, nside, cmask, nside/2)
                call udgrade_nest(tmp, nside, crhs, nside/2)
                
                ! W-stroke schedule
                do i = 1,2; call wstroke(nside/2, cmap, cmask, crhs); end do
                
                ! upgrade correction and update the map
                call udgrade_nest(cmap, nside/2, tmp, nside)
                forall (k=0:m) map(nn(1,k)) = map(nn(1,k)) + tmp(nn(1,k))
        end if
        
        ! post-smooth using Jacobi iteration
        do i = 1,niter
                forall (k=0:m) tmp(k) = (h2*rhs(nn(1,k)) - sum(L(2:9,k)*map(nn(2:9,k))))/L(1,k)
                forall (k=0:m) map(nn(1,k)) = tmp(k)
        end do
        
        ! calculate residual if requested
        if (present(residual)) then
                residual = 0.0
                forall (k=0:m) tmp(k) = rhs(nn(1,k)) - sum(L(:,k)*map(nn(:,k)))/h2
                forall (k=0:m) residual(nn(1,k)) = tmp(k)
                
                write (*,*) "Average residual", sqrt(sum(tmp(0:m)**2))/(m+1.0)
        end if
end subroutine wstroke

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! return nearest neighbours list (in arbitrary map ordering)
subroutine neighbours(nside, order, i, nn, k)
        integer nside, order, i, j, k, nn(8)
        
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
subroutine stencil(nside, order, i, nn, count, bartjan, laplace)
        integer nside, order, i, j, k, nn(9)
        integer, intent(out), optional :: count
        
        ! supported stencils - see below for description
        real(DP), dimension(9), intent(out), optional :: bartjan, laplace
        
        ! nearest neighbour list
        nn(1) = i; call neighbours(nside, order, i, nn(2:), k)
        if (k < 8) nn(9) = i; if (present(count)) count = k+1
        
        ! nearest neighbour average used by Bartjan van Tent et. al.
        if (present(bartjan)) then; bartjan = 0.0; bartjan(2:k+1) = 1.0/k; end if
        
        ! placeholder, this does not handle pixels with seven neighnours correctly
        if (present(laplace)) then
                laplace(1) = -8.0/3.0; laplace(2:9) = 1.0/3.0
                !if (k < 8) call warning("n=7 stencil not supported in Laplacian")
        end if
end subroutine stencil


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

end
