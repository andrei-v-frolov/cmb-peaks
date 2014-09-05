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

integer, parameter :: IO = SP                   ! default I/O precision
integer, parameter :: RING = 1, NEST = 2        ! ordering literals
logical, parameter :: verbose = .true.          ! diagnostic output

integer :: nmaps = 0, nside = 0, ord = 0, n = 0
character(len=80) :: header(64), fin1, op, fin2, fout
real(IO), dimension(:,:), allocatable :: M1, M2, Mout
logical, dimension(:,:), allocatable :: valid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type multigrid
        integer nside, n, m; real(DP) h2                ! grid specification
        real(DP), dimension(:), allocatable :: map, rhs, tmp
        real(DP), dimension(:,:), allocatable :: LAPL   ! Laplacian stencil
        integer,  dimension(:,:), allocatable :: nn     ! nearest neighbours
end type

type(multigrid), allocatable :: mg(:)

#define $MG(X) X => mg(l)%X
#define $MGVARS$ $MG(nside), $MG(n), $MG(m), $MG(h2), $MG(map), $MG(rhs), $MG(tmp), $MG(LAPL), $MG(nn)


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

! inpaint map using multigrid diffusion where mask is not unity
subroutine inpaint(map, mask, mout)
        real(IO), dimension(0:n) :: map, mask, mout
        integer i
        
        if (verbose) write (*,*) "Initalizing multigrid, masked pixel counts:"
        call mg_init(mg, nside, ord, map*mask, mask)
        
        associate(result => mg(1)%map, residual => mg(1)%tmp, m => mg(1)%m)
        if (verbose) write (*,*) "Running W-stroke iterations, average/max residual:"
        do i = 1,32
                call mg_wstroke(mg, 1)
                
                ! output residual
                if (.not. verbose) cycle
                call mg_residual(mg, 1, residual)
                write (*,*) sqrt(sum(residual**2))/m, maxval(abs(residual))
        end do
        
        mout = result
        end associate
        
        ! convert back to original order
        if (ord == RING) call convert_nest2ring(nside, mout)
end subroutine inpaint

! init multigrid structure
subroutine mg_init(mg, fside, order, imap, imask); use udgrade_nr
        type(multigrid), allocatable :: mg(:)
        integer i, k, l, fside, order, levels
        real(IO), dimension(0:12*fside**2-1) :: imap, imask
        
        ! total multigrid levels
        levels = log(nside/4.0)/log(2.0)
        if (levels < 1) levels = 1
        allocate(mg(levels))
        
        ! allocate grid storage
        do l = 1,levels; associate($MGVARS$)
                nside = ishft(fside,1-l)
                n = nside2npix(nside)-1
                h2 = (pi/3.0)/nside**2
                
                allocate(mg(l)%map(0:n), mg(l)%rhs(0:n), mg(l)%tmp(0:n), source=0.0)
                allocate(mg(l)%nn(9,0:n), mg(l)%LAPL(9,0:n))
        end associate; end do
        
        ! initialize finest grid
        l = 1; associate($MGVARS$, mask => mg(l)%tmp)
                map = imap; if (order == RING) call convert_ring2nest(nside, map)
                mask = imask; if (order == RING) call convert_ring2nest(nside, mask)
        end associate
        
        ! initialize coarse grids
        do l = 2,levels; associate($MGVARS$, mask => mg(l)%tmp)
                call udgrade_nest(mg(l-1)%tmp, mg(l-1)%nside, mask, nside)
        end associate; end do
        
        ! initalize stencils
        do l = 1,levels; associate($MGVARS$, mask => mg(l)%tmp)
                k = 0; do i = 0,n; if (mask(i) == 1.0) cycle
                        call stencil(nside, NEST, i, nn(:,k), Lw=LAPL(:,k)); k = k+1
                end do; m = k-1
                
                if (verbose) write (*,*) l, nside, m
        end associate; end do
end subroutine mg_init

! smooth map using Jacobi iteration
subroutine mg_smooth(mg, l, iterations)
        type(multigrid) mg(:); integer i, k, l, iterations
        
        associate($MGVARS$)
        do i = 1,iterations
                forall (k=0:m) tmp(k) = (h2*rhs(nn(1,k)) - sum(LAPL(2:9,k)*map(nn(2:9,k))))/LAPL(1,k)
                forall (k=0:m) map(nn(1,k)) = tmp(k)
        end do
        end associate
end subroutine mg_smooth

! calculate residual
subroutine mg_residual(mg, l, residual)
        type(multigrid) mg(:); integer i, k, l; real(DP) residual(:)
        
        residual = 0.0
        
        associate($MGVARS$)
        forall (k=0:m) residual(nn(1,k)) = rhs(nn(1,k)) - sum(LAPL(:,k)*map(nn(:,k)))/h2
        end associate
end subroutine mg_residual

! multigrid W-stroke: inpaints level-l map with L[map] = rhs where masked
recursive subroutine mg_wstroke(mg, l); use udgrade_nr
        type(multigrid) mg(:); integer i, k, l
        
        ! pre-smooth
        call mg_smooth(mg, l, 16)
        
        ! solve coarse problem
        if (l < size(mg)) then
                associate($MGVARS$, cmap => mg(l+1)%map, crhs => mg(l+1)%rhs)
                
                ! downgrade residual
                call mg_residual(mg, l, tmp)
                call udgrade_nest(tmp, nside, crhs, mg(l+1)%nside)
                
                ! W-stroke schedule
                cmap = 0.0; do i = 1,2; call mg_wstroke(mg, l+1); end do
                
                ! upgrade correction and update the map
                call udgrade_nest(cmap, mg(l+1)%nside, tmp, nside)
                forall (k=0:m) map(nn(1,k)) = map(nn(1,k)) + tmp(nn(1,k))
        end associate; else; call mg_smooth(mg, l, 96); end if
        
        ! post-smooth
        call mg_smooth(mg, l, 16)
end subroutine mg_wstroke


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
subroutine pix2gno(nside, order, i, p, XY, R)
        integer nside, order, i, p; real(DP) XY(2), R
        real(DP) theta, phi, U(3), V(3), W(3), X(3), Y(3)
        optional XY, R
        
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
        W = V/sum(U*V) - U; if (present(R)) R = sqrt(sum(W*W))
        
        ! bail unless we want coordinates
        if (.not. present(XY)) return
        
        ! project onto local orthonormal basis in (theta,phi) directions
        X = (/ cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta) /)
        Y = (/ -sin(phi), cos(phi), 0.0 /)
        
        XY = (/ sum(W*X), sum(W*Y) /)
end subroutine pix2gno

! calculate differential operator stencils (to be applied to nearest neighbours)
subroutine stencil(nside, order, i, nn, count, La, Lw, Lx, Dx, Dy, Dxx, Dxy, Dyy)
        integer nside, order, i, j, k, nn(9), count, status
        real(DP), dimension(9) :: La, Lw, Lx, Dx, Dy, Dxx, Dxy, Dyy
        
        ! interface description
        intent(in) nside, order, i
        intent(out) nn, count, La, Lw, Lx, Dx, Dy, Dxx, Dxy, Dyy
        optional count, La, Lw, Lx, Dx, Dy, Dxx, Dxy, Dyy
        
        ! working variables
        real, parameter :: gamma = -1.0
        integer, parameter :: m = 9, n = 10
        real(DP) x, y, XY(2), A(m), B(m), F(m,n), S(m), U(m,m), Q(n,m), V(n,n), W(m*n)
        
        ! nearest neighbour list
        nn(1) = i; call neighbours(nside, order, i, nn(2:), k)
        if (k < 8) nn(9) = i; if (present(count)) count = k+1
        
        ! average Laplacian stencil (used by Bartjan van Tent et. al.)
        if (present(La)) then
                W(1:k) = (8.0/3.0)/k; W(k+1:) = 0.0
                
                La = (/ -sum(W(1:k)), W(1:8) /)
        end if
        
        ! distance-weighted Laplacian stencil (better, still cheap to calculate)
        if (present(Lw)) then
                do j = 1,k; call pix2gno(nside, order, i, nn(j+1), R=S(j)); end do
                W(1:k) = 4.0*S(1:k)**gamma/sum(S(1:k)**(gamma+2.0)); W(k+1:) = 0.0
                
                Lw = (/ -sum(W(1:k)), W(1:8) /) * (pi/3.0)/nside**2
        end if
        
        ! bail unless exact (and expensive!) stencils are requested
        if (.not. any((/ present(Lx), present(Dx), present(Dy), present(Dxx), present(Dxy), present(Dyy) /))) return
        
        ! local polynomial basis (in gnomonic coordinates)
        do j = 1,m
                call pix2gno(nside, order, i, nn(j), XY)
                XY = nside/sqrt(pi/3.0) * XY; x = XY(1); y = XY(2)
                
                F(j,:) = (/ 1.0, x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y /)
        end do
        
        ! calculate singular value decomposition of basis matrix
        call dgesvd('All', 'All', m, n, F, m, S, U, m, V, n, W, size(W), status)
        if (status /= 0) call abort("SVD failed in stencil()")
        
        ! invert SVD to obtain stencil projector; last eigenvalue could be degenerate
        Q = 0.0; forall (j=1:k+1) Q(j,j) = 1.0/S(j)
        Q = matmul(Q,transpose(U))
        Q = matmul(transpose(V),Q)
        
        ! pixels 1 and 9 are the same in 7-neighbour case, roll 'em up together
        if (k < 8) then; Q(:,1) = Q(:,1) + Q(:,9); Q(:,9) = 0.0; end if
        
        ! return requested stencils
        if (present(Lx)) Lx = 2.0*(Q(4,:)+Q(6,:))
        if (present(Dx)) Dx = Q(2,:)
        if (present(Dy)) Dy = Q(3,:)
        if (present(Dxx)) Dxx = 2.0*Q(4,:)
        if (present(Dxy)) Dxy = Q(5,:)
        if (present(Dyy)) Dyy = 2.0*Q(6,:)
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