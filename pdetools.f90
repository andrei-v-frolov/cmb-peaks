module pdetools ! PDE operations on a HEALPix grid

use mapio
use udgrade_nr

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type grid
	integer nside, n, m; real(DP) h2                        ! grid specification
	real(DP), dimension(:), allocatable :: map, rhs, tmp    ! variables on a HEALPix grid (0:n)
	real(DP), dimension(:,:), allocatable :: laplacian      ! packed Laplacian stencil  (9,0:m)
	integer,  dimension(:,:), allocatable :: nn             ! packed nearest neighbours (9,0:m)
end type

#define $MG(X) X => mg(l)%X
#define $MGVARS$ $MG(nside), $MG(n), $MG(m), $MG(h2), $MG(map), $MG(rhs), $MG(tmp), $MG(laplacian), $MG(nn)

public :: inpaint, stencil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inpaint map using multigrid diffusion where mask is not unity
subroutine inpaint(nside, order, map, mask, mout, fill, apo)
	integer nside, order, i, mode
	real(IO), dimension(0:12*nside**2-1) :: map, mask, fill, apo, mout
	type(grid), allocatable :: mg(:); optional fill, apo
	
	if (verbose) write (*,*) "Initalizing multigrid, masked pixel counts:"
	
	! inpainting mode
	mode = 0
	if (present(fill)) mode = mode + 1
	if (present(apo))  mode = mode + 2
	
	select case (mode)
		case (3); call mg_init(mg, nside, order, map*mask, mask, fill, apo)
		case (2); call mg_init(mg, nside, order, map*mask, mask, map,  apo)
		case (1); call mg_init(mg, nside, order, map*mask, mask, fill)
		case (0); call mg_init(mg, nside, order, map*mask, mask)
	end select
	
	! run multigrid iterations
	associate(result => mg(1)%map, residual => mg(1)%tmp, m => mg(1)%m)
	if (verbose) write (*,*) "Running W-stroke iterations, average/max residual:"
	do i = 1,18
		! run W-stroke
		call mg_wstroke(mg, 1)
		
		! feather up on the last iteration
		if (i == 18) call mg_smooth(mg, 1, 8)
		
		! output residual
		if (.not. verbose) cycle
		call mg_residual(mg, 1)
		write (*,*) sqrt(sum(residual**2))/m, maxval(abs(residual))
	end do
	
	mout = result
	end associate
	
	! clean up
	call mg_free(mg)
	
	! convert back to original order
	if (order == RING) call convert_nest2ring(nside, mout)
end subroutine inpaint

! init multigrid structure
subroutine mg_init(mg, fside, order, imap, imask, fill, apo)
	type(grid), allocatable :: mg(:)
	integer i, k, l, fside, order, levels
	real(IO), dimension(0:12*fside**2-1) :: imap, imask, fill, apo
	optional fill, apo
	
	! total multigrid levels
	levels = log(fside/4.0)/log(2.0)
	if (levels < 1) levels = 1
	allocate(mg(levels))
	
	! allocate grid storage
	do l = 1,levels; associate($MGVARS$)
		nside = ishft(fside,1-l)
		n = nside2npix(nside)-1
		h2 = (pi/3.0)/nside**2
		
		allocate(mg(l)%map(0:n), mg(l)%rhs(0:n), mg(l)%tmp(0:n), source=0.0)
		allocate(mg(l)%nn(9,0:n), mg(l)%laplacian(9,0:n))
	end associate; end do
	
	! initialize finest grid
	l = 1; associate($MGVARS$)
		map = imap;  if (order == RING) call convert_ring2nest(nside, map)
		tmp = imask; if (order == RING) call convert_ring2nest(nside, tmp)
	end associate
	
	! initialize coarse grids
	do l = 2,levels; associate($MGVARS$)
		call udgrade_nest(mg(l-1)%tmp, mg(l-1)%nside, tmp, nside)
	end associate; end do
	
	! initialize stencils
	do l = 1,levels; associate($MGVARS$, mask => mg(l)%tmp)
		k = 0; do i = 0,n; if (mask(i) /= 0.0) cycle
			call stencil(nside, NEST, i, nn(:,k), Lw=laplacian(:,k)); k = k+1
		end do; m = k-1
		
		if (verbose) write (*,*) l, nside, m
	end associate; end do
	
	! initialize RHS with laplacian of a source map
	if (present(fill)) then; l = 1; associate($MGVARS$, src => mg(l)%tmp)
		src = fill;  if (order == RING) call convert_ring2nest(nside, src)
		forall (k=0:m) rhs(nn(1,k)) = sum(laplacian(:,k)*src(nn(:,k)))/h2
	end associate; end if
	
	! apodize RHS if apodization mask is provided
	if (present(apo)) then; l = 1; associate($MGVARS$, src => mg(l)%tmp)
		src = apo;   if (order == RING) call convert_ring2nest(nside, src)
		rhs = rhs*src
	end associate; end if
end subroutine mg_init

! deallocate multigrid structure
subroutine mg_free(mg)
	type(grid), allocatable :: mg(:); integer l
	
	! release multigrid storage
	do l = 1,size(mg)
		deallocate(mg(l)%map, mg(l)%rhs, mg(l)%tmp, mg(l)%laplacian, mg(l)%nn)
	end do
	
	! release multigrid structure
	deallocate(mg)
end subroutine mg_free

! smooth map using Jacobi iteration
subroutine mg_smooth(mg, l, iterations)
	type(grid) mg(:); integer i, k, l, iterations
	
	associate($MGVARS$)
	do i = 1,iterations
		forall (k=0:m) tmp(k) = (h2*rhs(nn(1,k)) - sum(laplacian(2:9,k)*map(nn(2:9,k))))/laplacian(1,k)
		forall (k=0:m) map(nn(1,k)) = tmp(k)
	end do
	end associate
end subroutine mg_smooth

! calculate residual
subroutine mg_residual(mg, l)
	type(grid) mg(:); integer k, l
	
	associate($MGVARS$); tmp = 0.0
	forall (k=0:m) tmp(nn(1,k)) = rhs(nn(1,k)) - sum(laplacian(:,k)*map(nn(:,k)))/h2
	end associate
end subroutine mg_residual

! multigrid W-stroke: inpaints level-l map with L[map] = rhs where masked
recursive subroutine mg_wstroke(mg, l)
	type(grid) mg(:); integer i, k, l
	
	! pre-smooth
	call mg_smooth(mg, l, 8)
	
	! solve coarse problem (with at least a few pixels)
	if (l < size(mg) .and. mg(l+1)%m > 16) then
		associate($MGVARS$, cmap => mg(l+1)%map, crhs => mg(l+1)%rhs)
		
		! downgrade residual
		call mg_residual(mg, l)
		call udgrade_nest(tmp, nside, crhs, mg(l+1)%nside)
		
		! W-stroke schedule
		cmap = 0.0; do i = 1,2; call mg_wstroke(mg, l+1); end do
		
		! upgrade correction and update the map
		call udgrade_nest(cmap, mg(l+1)%nside, tmp, nside)
		forall (k=0:m) map(nn(1,k)) = map(nn(1,k)) + tmp(nn(1,k))
	end associate; else; call mg_smooth(mg, l, max(mg(l)%m-15,48)); end if
	
	! post-smooth
	call mg_smooth(mg, l, 8)
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
subroutine pix2gno(nside, order, i, p, XY, R2)
	integer nside, order, i, p; real(DP) XY(2), R2
	real(DP) theta, phi, U(3), V(3), W(3), X(3), Y(3)
	optional XY, R2
	
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
	W = V/sum(U*V) - U; if (present(R2)) R2 = sum(W*W)
	
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
		do j = 1,k; call pix2gno(nside, order, i, nn(j+1), R2=S(j)); end do
		S(1:k) = S(1:k) * nside**2/(pi/3.0); W(1:k) = exp(-S(1:k)/1.61)
		W(1:k) = 4.0*W(1:k)/sum(S(1:k)*W(1:k)); W(k+1:) = 0.0
		
		Lw = (/ -sum(W(1:k)), W(1:8) /)
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

end module