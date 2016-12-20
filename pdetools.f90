module pdetools ! PDE operations on a HEALPix grid

use mapio
use almtools

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
	complex(DP), dimension(:,:), allocatable :: laplacian   ! packed Laplacian stencil  (9,0:m)
	integer,  dimension(:,:), allocatable :: nn             ! packed nearest neighbours (9,0:m)
end type

#define $MG(X) X => mg(l)%X
#define $MGVARS$ $MG(nside), $MG(n), $MG(m), $MG(h2), $MG(map), $MG(rhs), $MG(tmp), $MG(laplacian), $MG(nn)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! multigrid inpainter, real and complex variants, in single & double precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _dp, name ## _zs, name ## _zd; end interface

GENERIC(inpaint)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! wrapper routines for inpainting of common HEALPix CMB map formats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generic iterfaces are implemented using Fortran preprocessor
#define GENERIC(name) interface name; module procedure name ## _sp, name ## _dp; end interface

GENERIC(inpaint_qu)
GENERIC(inpaint_purified_qu)

public :: inpaint, inpaint_qu, inpaint_purified_qu, stencil

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! real maps
#define DATA real
#define GRID type(grid_dp)
#define STENCIL Lw

! single precision
#define XP SP
#define VARIANT(name) name ## _sp
#include 'multigrid.fin'
#include 'inpaint-qu.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _dp
#include 'multigrid.fin'
#include 'inpaint-qu.fin'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! complex maps
#define DATA complex
#define GRID type(grid_zd)
#define STENCIL Lz

! single precision
#define XP SP
#define VARIANT(name) name ## _zs
#include 'multigrid.fin'

! double precision
#define XP DP
#define VARIANT(name) name ## _zd
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
subroutine pix2gno(nside, order, i, p, XY, Z, R2, zeta)
	integer nside, order, i, p
	real(DP) XY(2), R2, zeta; complex(DPC) Z
	real(DP) theta, phi, alpha, beta, ua, ub
	real(DP) U(3), V(3), W(3), X(3), Y(3), A(3), B(3)
	optional XY, Z, R2, zeta
	
	! do we need to calculate stuff?
	logical needXY, needAB
	
	needAB = present(Z) .or. present(zeta)
	needXY = present(XY) .or. needAB
	
	! convert pixels to coordinates
	select case(order)
		case(RING)
			call pix2vec_ring(nside, i, U); if (needXY) call pix2ang_ring(nside, i, theta, phi)
			call pix2vec_ring(nside, p, V); if (needAB) call pix2ang_ring(nside, p, alpha, beta)
		case(NEST)
			call pix2vec_nest(nside, i, U); if (needXY) call pix2ang_nest(nside, i, theta, phi)
			call pix2vec_nest(nside, p, V); if (needAB) call pix2ang_nest(nside, p, alpha, beta)
		case default
			call abort(": ordering not supported")
	end select
	
	! project onto a tangent plane through origin
	W = V/sum(U*V) - U; if (present(R2)) R2 = sum(W*W)
	
	! bail unless we want coordinates
	if (.not. needXY) return
	
	! project onto local orthonormal basis in (theta,phi) directions
	X = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]
	Y = [-sin(phi), cos(phi), 0.0]
	
	if (present(XY)) XY = [sum(W*X), sum(W*Y)]
	
	! bail unless we want rotation
	if (.not. needAB) return
	
	! pixel p has different basis orientation, rotated by angle theta
	A = [cos(alpha)*cos(beta), cos(alpha)*sin(beta), -sin(alpha)]
	B = [-sin(beta), cos(beta), 0.0]
	
	! project A and B onto tangent plane
	!ua = sum(U*A); A = (A - ua*U)/sqrt(1.0 - ua*ua)
	!ub = sum(U*B); B = (B - ub*U)/sqrt(1.0 - ub*ub)
	
	! rotation angle is defined via 2x2 SVD
	theta = atan2(sum(A*Y-B*X),sum(A*X+B*Y))
	if (present(Z)) Z = exp((0.0,2.0)*theta)
	if (present(zeta)) zeta = theta
end subroutine pix2gno

! calculate differential operator stencils (to be applied to nearest neighbours)
subroutine stencil(nside, order, i, nn, count, La, Lw, Lz, Lx, Dx, Dy, Dxx, Dxy, Dyy)
	integer nside, order, i, j, k, nn(9), count, status
	real(DP), dimension(9) :: La, Lw, Lx, Dx, Dy, Dxx, Dxy, Dyy
	complex(DPC), dimension(9) :: Lz
	
	! interface description
	intent(in) nside, order, i
	intent(out) nn, count, La, Lw, Lz, Lx, Dx, Dy, Dxx, Dxy, Dyy
	optional count, La, Lw, Lz, Lx, Dx, Dy, Dxx, Dxy, Dyy
	
	! working variables
	integer, parameter :: m = 9, n = 10
	real(DP) x, y, XY(2), A(m), B(m), F(m,n), S(m), U(m,m), Q(n,m), V(n,n), W(m*n)
	complex(DPC) Z(8)
	
	! nearest neighbour list
	nn(1) = i; call neighbours(nside, order, i, nn(2:), k)
	if (k < 8) nn(9) = i; if (present(count)) count = k+1
	
	! average Laplacian stencil (used by Bartjan van Tent et. al.)
	if (present(La)) then
		W(1:k) = (8.0/3.0)/k; W(k+1:) = 0.0
		
		La = [-sum(W(1:k)), W(1:8)]
	end if
	
	! distance-weighted Laplacian stencil (better, still cheap to calculate)
	if (present(Lw)) then
		do j = 1,k; call pix2gno(nside, order, i, nn(j+1), R2=S(j)); end do
		S(1:k) = S(1:k) * nside**2/(pi/3.0); W(1:k) = exp(-S(1:k)/1.61)
		W(1:k) = 4.0*W(1:k)/sum(S(1:k)*W(1:k)); W(k+1:) = 0.0
		
		Lw = [-sum(W(1:k)), W(1:8)]
	end if
	
	! distance-weighted complex Laplacian stencil for QU maps
	if (present(Lz)) then
		do j = 1,k; call pix2gno(nside, order, i, nn(j+1), R2=S(j), Z=Z(j)); end do
		S(1:k) = S(1:k) * nside**2/(pi/3.0); W(1:k) = exp(-S(1:k)/1.61)
		W(1:k) = 4.0*W(1:k)/sum(S(1:k)*W(1:k)); W(k+1:) = 0.0
		
		Lz = [cmplx(-sum(W(1:k))), Z * W(1:8)]
	end if
	
	! bail unless exact (and expensive!) stencils are requested
	if (.not. any([present(Lx), present(Dx), present(Dy), present(Dxx), present(Dxy), present(Dyy)])) return
	
	! local polynomial basis (in gnomonic coordinates)
	do j = 1,m
		call pix2gno(nside, order, i, nn(j), XY)
		XY = nside/sqrt(pi/3.0) * XY; x = XY(1); y = XY(2)
		
		F(j,:) = [1.0, x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y]
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