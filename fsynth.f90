! $Id$
! Synthesize a filtered map, applying a mask correction, and output hot pixels
! invoke: wiener ... | fsynth <mmap-alms.fits> <map.fits[:iqu]> [mask-alms.fits[:tdb]]

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

integer, parameter :: lmax = 1535, nside = 512, hmax = 256
character(len=80) :: header(hmax), fmap, fmask, fout
real(DP), allocatable :: Mmap(:,:), Mout(:,:), Mask(:,:)
real(DP) beam(0:lmax)

integer i, l, n, nmaps, mmaps, io; real :: b, tol = 1.0

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
call getArgument(2, fout)

i = index(fout, ":iqu", .true.); if (i > 0) fout(i:) = ""
call make_map(fmap, beam, lmax, nside, Mmap, n, nmaps, (i > 0))

! make filtered mask
if (nArguments() < 3) then
	allocate (Mask(0:n,1)); Mask = 1.0; mmaps = 1
else
	call getArgument(3, fmask); i = index(fmask, ":", .true.)
	if (i > 0) then; read (fmask(i+1:),*) tol; fmask(i:) = ""; end if
	
	call make_map(fmask, beam/beam(0), lmax, nside, Mask, n, mmaps, .false.)
	
	! remove non-compliant pixels
	where (Mask <= 0.0 .or. abs(log10(Mask)) > tol/10.0) Mask = 0.0
end if

! mask correction
allocate (Mout(0:n,nmaps))
do i = 0,n
	Mout(i,:) = Mmap(i,:)/Mask(i,1)
end do

! output extrema distribution in first extent
call extrema(Mout(:,1), Mask(:,1), nside, n)

! output corrected map
call write_minimal_header(header, 'MAP', nside=nside, order=2, creator='FSYNTH', version='$Revision$', polar=(nmaps .eq. 3))
call output_map(Mout, header, '!'//fout)


contains

! synthesize a map from alms convolved with a beam function
subroutine make_map(falms, bls, lmax, nside, map, n, nout, decompose)
	character(*) falms; real bls(0:lmax)
	integer lmax, nalms, nside, next, nout, j, k, l, m, n
	logical, optional :: decompose; logical iqu
	character(len=80) :: header(hmax,3)
	real(DP), allocatable :: aks(:,:,:), map(:,:)
	complex(DPC), allocatable :: alms(:,:,:)
	
	! read in alms data
	iqu = present(decompose) .and. decompose
	nalms = number_of_alms(falms, next); n = nside2npix(nside)-1
	if (iqu) then; next = 1; nout = 3; else; nout = next; end if
	allocate(aks(nalms,4,next), alms(next,0:lmax,0:lmax), map(0:n,nout))
	
	call fits2alms(falms, nalms, aks, 3, header, hmax, next)
	
	! convert to complex alms, applying beam
	do j = 1,next; do k = 1,nalms
		l = aks(k,1,j); m = aks(k,2,j)
		if (l > lmax .or. m > lmax) cycle
		alms(j,l,m) = bls(l)*cmplx(aks(k,3,j), aks(k,4,j))
	end do; end do
	
	! synthesize map (in nested ordering)
	if (iqu) then
		call alm2map_iqu(nside, lmax, lmax, alms, map(:,1:nout))
	else
		if (next == 1) call alm2map(nside, lmax, lmax, alms, map(:,1))
		if (next >  1) call alm2map(nside, lmax, lmax, alms, map(:,1:next))
	end if
	call convert_ring2nest(nside, map)
	
	deallocate(aks, alms)
end subroutine

! synthesize an IQU decomposition map of tidal tensor (using canned HEALPix derivatives)
subroutine alm2map_iqu_hpx(nside, lmax, mmax, alms, Miqu)
	integer p, l, m, n, nside, lmax, mmax
	complex(DPC) alms(1,0:lmax,0:mmax)
	real(DP) Miqu(12*nside**2,3), theta, phi
	
	! temporary storage
	real(DP), allocatable :: I(:), C(:), D1(:,:), D2(:,:)
	complex(DPC), allocatable :: tlms(:,:,:)
	
	n = nside2npix(nside)-1
	allocate(I(0:n), C(0:n), D1(0:n,2), D2(0:n,3), tlms(1,0:lmax,0:mmax))
	
	! compute alm's of derivative maps
	do l = 1,lmax; do m = 0,min(l,mmax)
		tlms(1,l,m) = -alms(1,l,m)/(l+1.0)/l
	end do; end do
	 
	! synthesize component maps from alm's
	call alm2map(nside, lmax, mmax, alms, I)
	call alm2map_der(nside, lmax, mmax, tlms, C, D1, D2)

	! synthesize cot(theta) map (in ring ordering)
	do p = 0,n
		call pix2ang_ring(nside, p, theta, phi); C(p) = cotan(theta)
	end do
	
	! assemble output map from components
	Miqu(:,1) = I
	Miqu(:,2) = D2(:,1) - D2(:,3) - C*D1(:,1)
	Miqu(:,3) = 2.0*(D2(:,2) - C*D1(:,2))
	
	! clean up after ourselves
	deallocate(tlms, I, C, D1, D2)
end subroutine

! synthesize an IQU decomposition map of tidal tensor (using a_lm coefficient recursion)
subroutine alm2map_iqu(nside, lmax, mmax, alms, Miqu)
	integer p, l, m, n, nside, lmax, mmax
	complex(DPC) alms(1,0:lmax,0:mmax), mu
	real(DP) Miqu(12*nside**2,3), theta, phi
	
	! temporary storage
	real(DP), allocatable :: I(:), T(:), Q(:), U(:), S(:)
	complex(DPC), allocatable :: tlms(:,:,:), qlms(:,:,:), ulms(:,:,:)
	
	n = nside2npix(nside)-1
	allocate(I(0:n), T(0:n), Q(0:n), U(0:n), S(0:n))
	allocate(tlms(1,0:lmax,0:mmax), qlms(1,0:lmax,0:mmax), ulms(1,0:lmax,0:mmax))
	
	tlms = 0.0; qlms = 0.0; ulms = 0.0
	
	! compute alm's of derivative maps
	do l = 2,lmax; do m = 0,min(l,mmax)
		mu = (0.0,2.0)/(l-1.0) * m/l
		
		tlms(1,l,m) = (l-1.0)/(l+1.0) * alms(1,l,m)
		qlms(1,l,m) = (l-2*m*m)/(l-0.5)/l * tlms(1,l,m)
		ulms(1,l,m) = -mu*(l-2) * A(l,m) * alms(1,l-1,m)
		if (m < l-1) ulms(1,l-2,m) = ulms(1,l-2,m) + mu*(l+1) * A(l-1,m) * alms(1,l-1,m)
		if (m < l-1) qlms(1,l-2,m) = qlms(1,l-2,m) - 4.0/l*(l+0.5)/(l+1.0) * A(l,m)*A(l-1,m) * alms(1,l,m)
	end do; end do
	
	! synthesize component maps from alm's
	call alm2map(nside, lmax, mmax, alms, I)
	call alm2map(nside, lmax, mmax, tlms, T)
	call alm2map(nside, lmax, mmax, qlms, Q)
	call alm2map(nside, lmax, mmax, ulms, U)
	
	! synthesize sin(theta)^2 map (in ring ordering)
	do p = 0,n
		call pix2ang_ring(nside, p, theta, phi); S(p) = sin(theta)**2
	end do
	
	! assemble output map from components
	Miqu(:,1) = I; Miqu(:,2) = Q/S + T; Miqu(:,3) = U/S
	
	! clean up after ourselves
	deallocate(tlms, qlms, ulms, I, T, Q, U, S)
end subroutine

! shorthand for normalization pre-factor
function A(l,m)
	integer l, m; real A
	
	A = sqrt((l*l-m*m)/(4.0*l*l-1.0))
end function

! find extrema in the interior of the masked map
subroutine extrema(M, mask, nside, n)
	integer i, k, n, nside, nn(8)
	real(DP), dimension(0:n) :: M, mask
	
	do i = 0,n
		if (mask(i) == 0.0) cycle
		call neighbours_nest(nside, i, nn, k)
		if (any(mask(nn(1:k)) == 0.0)) cycle
		
		if (M(i) > maxval(M(nn(1:k)))) write (*,'(I,G24.16,I)') i, M(i), +1
		if (M(i) < minval(M(nn(1:k)))) write (*,'(I,G24.16,I)') i, M(i), -1
	end do
end subroutine

end
