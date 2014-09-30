! $Id$
! create a blank map from pixel listing using existing map format
! pixel values are to be supplied to stdin in triples 'theta phi value'
! invoke: pxl2map <map.fits> <out.fits>

program pxl2map

! HEALPix includes
use mapio
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: fin, fout
real(IO), allocatable :: M(:,:)
integer, allocatable, target :: idx(:)
integer :: nside = 0, nmaps = 0, ord = 0, n, status

real(DP) theta, phi, value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! template and output maps
call getArgument(1, fin )
call getArgument(2, fout)

! import input map
call read_map(fin, M, nside, nmaps, ord)
n = nside2npix(nside)-1; allocate(idx(0:n))

M = 1.0/0.0

! input loop
do
	read (*,*,iostat=status) theta, phi, value
	if (status < 0) exit
	M(pixel(theta,phi),:) = value
	!M(disk(theta,phi,90.0),:) = value
	!M(disk(theta,phi,600.0),:) = M(disk(theta,phi,600.0),:) + 1.0
end do

! output map
call write_map(fout, M, nside, ord, creator='PXL2MAP')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nearest pixel index
function pixel(theta, phi)
	real(DP) theta, phi, v(3); integer pixel
	
	select case (ord)
		case(RING); call ang2pix_ring(nside, theta, phi, pixel)
		case(NEST); call ang2pix_nest(nside, theta, phi, pixel)
		case default; call abort(": ordering not supported")
	end select
end function

! disk pixel index 
function disk(theta, phi, radius)
	real(DP) theta, phi, radius, v(3)
	integer, pointer :: disk(:); integer k
	
	call ang2vec(theta, phi, v)
	call query_disc(nside, v, radius * pi/(180.0*60.0), idx, k, nest=ord-1)
	
	disk => idx(0:k-1)
end function

end
