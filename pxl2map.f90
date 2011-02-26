! $Id$
! create a blank map from pixel listing using existing map format
! pixel values are to be supplied to stdin in pairs 'pixel value'
! invoke: pxl2map <map.fits> <out.fits>

program pxl2map

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: header(64), fin, fout
integer nmaps, nside, npix, n, ntot, ord, io
real(SP), allocatable :: M(:,:)

integer i; real(DP) p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! template and output maps
call getArgument(1, fin )
call getArgument(2, fout)

! import input map
ntot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
npix = nside2npix(nside); n = npix-1

allocate(M(0:n,nmaps))

call input_map(fin, M, npix, nmaps)

M = 0.0

! input loop
do
	read (*,*,iostat=io) i, p
	if (io < 0) exit
	M(i,:) = p
end do

call write_minimal_header(header, 'MAP', nside=nside, order=ord)
call output_map(M, header, '!'//fout)

end
