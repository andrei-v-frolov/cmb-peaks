! $Id$
! Synthesize a filtered map, applying a mask correction
! invoke: wiener ... | fsynth <mmap-alms.fits> <map.fits> <mask-alms.fits[:tdb]>

program fsynth

! HEALPix includes
use extension
use head_fits
use healpix_types
use fitstools
use alm_tools
use pix_tools

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: lmax = 1024, nside = 512
character(len=80) :: header(64), fmap, fmask, fout
real(DP), allocatable :: Mmap(:,:), Mout(:,:), Mask(:,:)
real(DP) beam(0:lmax)

integer i, l, n, nmaps, io; real :: b, tol = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! read beam from stdin
beam = 0.0; do
	read (*,*,iostat=io) l, b
	if (io < 0) exit
	if (l > lmax) cycle
	beam(l) = b
end do

! make filtered map
call getArgument(1, fmap)
call make_map(fmap, beam, lmax, nside, Mmap, n, nmaps)

! make filtered mask
if (nArguments() < 3) then
	allocate (Mask(0:n,nmaps)); Mask = 1.0
else
	call getArgument(3, fmask); i = index(fmask, ":", .true.)
	if (i > 0) then; read (fmask(i+1:),*) tol; fmask(i:) = ""; end if
	
	call make_map(fmask, beam/beam(0), lmax, nside, Mask, n, nmaps)
	
	! remove non-compliant pixels
	where (Mask <= 0.0 .or. abs(log10(Mask)) > tol/10.0) Mask = 0.0
end if

! mask correction
allocate (Mout(0:n,nmaps))
Mout = Mmap/Mask

! output corrected map
call getArgument(2, fout)
call write_minimal_header(header, 'MAP', nside=nside, order=1, creator='FSYNTH', version='$Revision$')
call output_map(Mout, header, '!'//fout)


contains

! synthesize a map from alms convolved with a beam function
subroutine make_map(falms, bls, lmax, nside, map, n, next)
	character(*) falms; real bls(0:lmax)
	integer lmax, nalms, nside, next, k, l, m, n
	character(len=80) :: header(64,3)
	real(DP), allocatable :: aks(:,:,:), map(:,:)
	complex(DPC), allocatable :: alms(:,:,:)
	
	nalms = number_of_alms(falms, next); n = nside2npix(nside)-1
	allocate(aks(nalms,4,next), alms(next,0:lmax,0:lmax), map(0:n,next))
	
	call fits2alms(falms, nalms, aks, 3, header, 64, next)
	
	do k = 1,nalms
		l = aks(k,1,1); m = aks(k,2,1)
		if (l > lmax .or. m > lmax) cycle
		alms(:,l,m) = bls(l)*cmplx(aks(k,3,:), aks(k,4,:))
	end do
	
	call alm2map(nside, lmax, lmax, alms, map(:,next))
	
	deallocate(aks, alms)
end subroutine

end
