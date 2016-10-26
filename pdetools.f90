module pdetools ! PDE operations on a HEALPix grid

use mapio
use udgrade_nr
use complex_qu

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type grid_dp
	integer nside, n, m; real(DP) h2                        ! grid specification
	real(DP), dimension(:), allocatable :: map, rhs, tmp    ! variables on a HEALPix grid (0:n)
	real(DP), dimension(:,:), allocatable :: laplacian      ! packed Laplacian stencil  (9,0:m)
	integer,  dimension(:,:), allocatable :: nn             ! packed nearest neighbours (9,0:m)
end type

type grid_zd
	integer nside, n, m; real(DP) h2                        ! grid specification
	complex(DP), dimension(:), allocatable :: map, rhs, tmp ! variables on a HEALPix grid (0:n)
	real(DP), dimension(:,:), allocatable :: laplacian      ! packed Laplacian stencil  (9,0:m)
	integer,  dimension(:,:), allocatable :: nn             ! packed nearest neighbours (9,0:m)
end type

#define $MG(X) X => mg(l)%X
#define $MGVARS$ $MG(nside), $MG(n), $MG(m), $MG(h2), $MG(map), $MG(rhs), $MG(tmp), $MG(laplacian), $MG(nn)

public :: inpaint, stencil

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _zs; end interface

GENERIC(inpaint)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! single precision, real maps
#define XP SP
#define DATA real
#define GRID type(grid_dp)
#define VARIANT(name) name ## _sp
#include 'multigrid.fin'

! single precision, complex maps
#define XP SP
#define DATA complex
#define GRID type(grid_zd)
#define VARIANT(name) name ## _zs
#include 'multigrid.fin'

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