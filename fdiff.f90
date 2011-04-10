! $Id$

program fdiff

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: header(64), fin1, fin2, fout
integer nmaps, nside, npix, ntot, ord
real(SP), allocatable :: M1(:,:), M2(:,:), Mout(:,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! input map
call getArgument(1, fin1 )
call getArgument(2, fin2 )
call getArgument(3, fout )

! import input map
ntot = getsize_fits(fin1, nmaps=nmaps, nside=nside, ordering=ord)

npix = nside2npix(nside)

allocate(M1(npix,nmaps), M2(npix,nmaps), Mout(npix,nmaps))

call input_map(fin1, M1, npix, nmaps)
call input_map(fin2, M2, npix, nmaps)

Mout = M1/M2

call write_minimal_header(header, 'MAP', nside=nside, order=ord)
call output_map(Mout, header, '!'//fout)

end
