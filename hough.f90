! $Id$
! HEALPix Hough transform, produces distribution of arcs tangent to supplied normal vectors
! hough V.fits nside output.fits

program vhist

! HEALPix includes
use mapio
use pix_tools
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=8000) :: fin, bins, fout
integer :: nmaps = 0, nside = 0, ord = 0, vec = HLPX
real(IO), allocatable :: M(:,:), Mout(:,:)
real(DP), allocatable :: hist(:,:)

integer i, hside, layers, status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
call getArgument(1, fin); call read_map(fin, M, nside, nmaps, ord, vec=vec)
if (nmaps < 2 .or. nmaps > 3) call abort("normal vector map is required as input")

call getArgument(2, bins); read (bins,*,iostat=status) hside
if (status /= 0) call abort("cannot parse nside in " // trim(bins))

call getArgument(3, fout)

! radial layer spacing should be quite fine
layers = min(256,4*hside)

! initialize accumulators
allocate(hist(0:12*hside**2-1,layers), Mout(0:12*hside**2-1,2)); hist = 0.0

! accumulate histogram
do i = 1,12*nside**2-1
	call accumulate(nside, ord, i, real(M(i,1:2),DP))
end do

! reduce the results
Mout(:,1) = maxval(hist,2)
Mout(:,2) = maxloc(hist,2)

! write output histogram
call write_map(fout, Mout, hside, RING, creator='HOUGH')

contains

subroutine accumulate(nside, order, p, normal)
	integer nside, order, p, q, i
	real(DP) alpha, w, normal(2), U(3), V(3), A(3), B(3)
	intent(in) nside, order, p, normal
	
	! convert pixel index to Cartesian coordinates
	select case(order)
		case(RING); call pix2vec_ring(nside, p, U)
		case(NEST); call pix2vec_nest(nside, p, U)
		case default; call abort(": ordering not supported")
	end select
	
	! convert normal vector to Cartesian coordinates
	V = hlpx2cart(nside, order, p, [normal/norm2(normal), 0.0])
	
	! scan possible intercept values
	do i = 1,layers
		alpha = (i-0.5)/layers * (pi/2)
		A = sin(alpha)*U; B = cos(alpha)*V; w = sum(normal**2)
		call vec2pix_ring(hside, A+B, q); hist(q,i) = hist(q,i) + w
		call vec2pix_ring(hside, A-B, q); hist(q,i) = hist(q,i) + w
	end do
end subroutine

end
