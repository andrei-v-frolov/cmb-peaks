! $Id$
! create a blank map from pixel listing using existing map format
! pixel values are to be supplied to stdin in triples 'theta phi value'
! invoke: pxl2map <map.fits> <out.fits>

program pxl2map

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: IO = SP               ! default I/O precision
integer, parameter :: RING = 1, NEST = 2    ! ordering literals

character(len=80) :: header(64), fin, fout
integer nmaps, nside, npix, n, ntot, ord, status
real(IO), allocatable :: M(:,:)
integer, allocatable, target :: idx(:)

real(DP) theta, phi, value


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! template and output maps
call getArgument(1, fin )
call getArgument(2, fout)

! import input map
ntot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
npix = nside2npix(nside); n = npix-1

allocate(M(0:n,nmaps), idx(npix))

!call input_map(fin, M, npix, nmaps)

!M = 0.0
M = 1.0/0.0

! input loop
do
	read (*,*,iostat=status) theta, phi, value
	if (status < 0) exit
	M(pixel(theta,phi),:) = value
	!M(disk(theta,phi,90.0),:) = value
	!M(disk(theta,phi,600.0),:) = M(disk(theta,phi,600.0),:) + 1.0
end do

call write_minimal_header(header, 'MAP', nside=nside, order=ord)
call output_map(M, header, '!'//fout)

contains

! nearest pixel index
function pixel(theta, phi)
	real(DP) theta, phi, warcmin, v(3); integer pixel
	
	select case (ord)
		case(RING); call ang2pix_ring(nside, theta, phi, pixel)
		case(NEST); call ang2pix_nest(nside, theta, phi, pixel)
		case default; call abort(": ordering not supported")
	end select
	
	write (*,*) pixel
end function

! disk pixel index 
function disk(theta, phi, warcmin)
	real(DP) theta, phi, warcmin, v(3)
	integer, pointer :: disk(:); integer k
	
	call ang2vec(theta, phi, v)
	call query_disc(nside, v, warcmin * pi/(180.0*60.0), idx, k, nest=ord-1)
	
	disk => idx(1:k)
end function

end
