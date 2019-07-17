! $Id$
! HEALPix vector histogram generator, produces distribution of 3-vector directions in the input map
! vhist M.fits nside output.fits

program vhist

! HEALPix includes
use mapio
use pix_tools
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=8000) :: fin, bins, fout
integer :: nmaps = 0, nside = 0, ord = 0, pol = -1
real(IO), allocatable :: M(:,:), hist(:,:)

integer i, p, hside, status
real(DP) v(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
call getArgument(1, fin); call read_map(fin, M, nside, nmaps, ord, pol)
if (nmaps /= 3) call abort("3D vector map (in cartesian coordinates) is required as input")

call getArgument(2, bins); read (bins,*,iostat=status) hside
if (status /= 0) call abort("cannot parse nside in " // trim(bins))

call getArgument(3, fout)

! initialize accumulators
allocate(hist(0:12*hside**2-1,4)); hist = 0.0

! accumulate histogram
do i = 1,12*nside**2-1
	v = M(i,:)/norm2(M(i,:))
	call vec2pix_nest(hside, v, p)
	
	hist(p,1) = hist(p,1) + 1.0
	hist(p,2:4) = hist(p,2:4) + M(i,:)
end do

! write output histogram
call write_map(fout, hist, hside, NEST, creator='VHIST')

end
