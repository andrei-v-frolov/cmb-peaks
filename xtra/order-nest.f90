! $Id$

program ordernest

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none

integer :: nside = 64, ord = 2
character(len=80) :: header(64)
real(SP), allocatable :: M(:,:)
integer i, npix
real(DP) theta,phi

npix = nside2npix(nside)

allocate(M(npix,1))

do i = 1,npix
        M(i,1) = (i-1.0)/(npix-1.0)
        call pix2ang_nest(nside, i-1, theta, phi)
        write (*,*) theta,phi, i-1
end do

! write output map
call write_minimal_header(header, 'MAP', nside=nside, order=ord)
call output_map(M, header, '!order-nest.fits')

end
