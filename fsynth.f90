! $Id$
! Synthesize a filtered map, applying a mask correction, and output hot pixels
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

! output extrema distribution in first extent
call extrema(Mout(:,1), Mask(:,1), nside, n)

! output corrected map
call getArgument(2, fout)
call write_minimal_header(header, 'MAP', nside=nside, order=2, creator='FSYNTH', version='$Revision$')
call output_map(Mout, header, '!'//fout)


contains

! synthesize a map from alms convolved with a beam function
subroutine make_map(falms, bls, lmax, nside, map, n, next)
	character(*) falms; real bls(0:lmax)
	integer lmax, nalms, nside, next, j, k, l, m, n
	character(len=80) :: header(64,3)
	real(DP), allocatable :: aks(:,:,:), map(:,:)
	complex(DPC), allocatable :: alms(:,:,:)
	
	! read in alms data
	nalms = number_of_alms(falms, next); n = nside2npix(nside)-1
	allocate(aks(nalms,4,next), alms(next,0:lmax,0:lmax), map(0:n,next))
	
	call fits2alms(falms, nalms, aks, 3, header, 64, next)
	
	! convert to complex alms, applying beam
	do j = 1,next; do k = 1,nalms
		l = aks(k,1,j); m = aks(k,2,j)
		if (l > lmax .or. m > lmax) cycle
		alms(j,l,m) = bls(l)*cmplx(aks(k,3,j), aks(k,4,j))
	end do; end do
	
	! synthesize map (in nested ordering)
	if (next == 1) call alm2map(nside, lmax, lmax, alms, map(:,1))
	if (next >  1) call alm2map(nside, lmax, lmax, alms, map(:,1:next))
	call convert_ring2nest(nside, map)
	
	deallocate(aks, alms)
end subroutine

! find extrema in the interior of the masked map
subroutine extrema(M, mask, nside, n)
	integer i, k, n, nside, nn(8)
	real(DP), dimension(0:n) :: M, mask
	
	do i = 0,n
		if (mask(i) == 0.0) cycle
		call neighbours_nest(nside, i, nn, k)
		if (any(mask(nn(1:k)) == 0.0)) cycle
		
		if (M(i) > maxval(M(nn(1:k))) .or. M(i) < minval(M(nn(1:k)))) write (*,'(I,G24.16)') i, M(i)
	end do
end subroutine

end
