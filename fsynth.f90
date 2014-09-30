! $Id$
! Synthesize a filtered map, applying a mask correction, and output map statistics
! invoke: wiener ... | fsynth <mmap-alms.fits> <map.fits[:mode]> [mask-alms.fits[:tdb]]
!   - mode :iqu synthesizes IQU decomposition map of tidal tensor, outputs hot pixels
!   - mode :inv synthesizes curvature invariants map, outputs Minkowski functionals

program fsynth

! HEALPix includes
use extension
use alm_tools
use mapio

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: default = 0, iqu = 1, inv = 2; integer :: mode = default

character(len=80) :: fmap, fmask, fout
integer, parameter :: lmax = 4000, nside = 2048, hmax = 256
real(DP), allocatable :: mmap(:,:), mask(:,:), map(:,:)
real(IO), allocatable :: mout(:,:)
real(DP) beam(0:lmax)

integer i, l, n, nmaps, mmaps, status; real :: b, tol = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! read beam from stdin
beam = 0.0; do
	read (*,*,iostat=status) l, b
	if (status < 0) exit
	if (l > lmax) cycle
	beam(l) = b
end do

! make filtered map
call getArgument(1, fmap)
call getArgument(2, fout)

i = index(fout, ":iqu", .true.); if (i > 0) then; fout(i:) = ""; mode = iqu; end if
i = index(fout, ":inv", .true.); if (i > 0) then; fout(i:) = ""; mode = inv; end if
call make_map(fmap, beam, lmax, nside, mmap, n, nmaps, mode)

! make filtered mask
if (nArguments() < 3) then
	allocate (mask(0:n,1)); mask = 1.0; mmaps = 1
else
	call getArgument(3, fmask); i = index(fmask, ":", .true.)
	if (i > 0) then; read (fmask(i+1:),*) tol; fmask(i:) = ""; end if
	
	call make_map(fmask, beam/beam(0), lmax, nside, mask, n, mmaps, default)
	
	! remove non-compliant pixels
	where (mask <= 0.0 .or. abs(log10(mask)) > tol/10.0) mask = 0.0
end if

! mask correction
allocate (map(0:n,nmaps)); forall (i = 0:n) map(i,:) = mmap(i,:)/mask(i,1)

! output map statistics
select case (mode)
	case (default,iqu); call extrema(map(:,1), mask(:,1), nside, n)	! extrema distribution
	case (inv);         call minkowski(map, mask(:,1), nmaps, n, 1024)	! cumulative Minkowski functionals
	!case (inv);         call skeleton(map(:,4), mask(:,1), nside, n)	! ...
end select

! output corrected map (in specified precision)
allocate (mout(0:n,nmaps)); mout = map
call write_map(fout, mout, nside, NEST, creator='FSYNTH')

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! synthesize a map from alms convolved with a beam function
subroutine make_map(falms, bls, lmax, nside, map, n, nout, mode)
	character(*) falms; real bls(0:lmax)
	integer lmax, nalms, nside, next, nout, mode, j, k, l, m, n
	character(len=80) :: header(hmax,3)
	real(DP), allocatable :: aks(:,:,:), map(:,:)
	complex(DPC), allocatable :: alms(:,:,:)
	
	! read in alms data
	nalms = number_of_alms(falms, next); n = nside2npix(nside)-1
	select case (mode)
		case (iqu,inv); next = 1; nout = 3
		case default; nout = next
	end select
	allocate(aks(nalms,4,next), alms(next,0:lmax,0:lmax), map(0:n,nout))
	
	call fits2alms(falms, nalms, aks, 3, header, hmax, next)
	
	! convert to complex alms, applying beam
	do j = 1,next; do k = 1,nalms
		l = aks(k,1,j); m = aks(k,2,j)
		if (l > lmax .or. m > lmax) cycle
		alms(j,l,m) = bls(l)*cmplx(aks(k,3,j), aks(k,4,j))
	end do; end do
	
	! synthesize map (in nested ordering)
	select case (mode)
		case (iqu); call alm2map_iqu(nside, lmax, lmax, alms, map(:,1:nout))
		case (inv); call alm2map_inv(nside, lmax, lmax, alms, map(:,1:nout))
		case default;
			if (next == 1) call alm2map(nside, lmax, lmax, alms, map(:,1))
			if (next >  1) call alm2map(nside, lmax, lmax, alms, map(:,1:next))
	end select
	call convert_ring2nest(nside, map)
	
	deallocate(aks, alms)
end subroutine

! synthesize curvature invariants map for Minkowski functional evaluation (using canned HEALPix derivatives)
subroutine alm2map_inv(nside, lmax, mmax, alms, Minv)
	integer p, l, m, n, nside, lmax, mmax
	complex(DPC) alms(1,0:lmax,0:mmax)
	real(DP) Minv(0:12*nside**2-1,3), theta, phi
	
	! temporary storage
	real(DP), allocatable :: I(:), C(:), D1(:,:), D2(:,:)
	
	n = nside2npix(nside)-1
	allocate(I(0:n), C(0:n), D1(0:n,2), D2(0:n,3))
	
	! synthesize component maps from alm's
	call alm2map_der(nside, lmax, mmax, alms, I, D1, D2)
	
	do p = 0,n
		! cot(theta) map (in ring ordering)
		call pix2ang_ring(nside, p, theta, phi); C(p) = cotan(theta)
		
		! second covariant dervatives
		D2(p,2) = D2(p,2) - C(p) * D1(p,2)
		D2(p,3) = D2(p,3) + C(p) * D1(p,1)
		
		! covariant gradient squared
		C(p) = D1(p,1)**2 + D1(p,2)**2
		
		! Minkowski functional maps
		Minv(p,1) = I(p)
		Minv(p,2) = sqrt(C(p))
		Minv(p,3) = - (D2(p,1)*D1(p,2)**2 - 2.0*D2(p,2)*D1(p,1)*D1(p,2) + D2(p,3)*D1(p,1)**2)/C(p)
		
		! skeleton map
		!Minv(p,4) = (D2(p,3) - D2(p,1))*D1(p,1)*D1(p,2) + D2(p,2) * (D1(p,1)**2 - D1(p,2)**2)
	end do
	
	! clean up after ourselves
	deallocate(I, C, D1, D2)
end subroutine

! synthesize an IQU decomposition map of tidal tensor (using canned HEALPix derivatives)
subroutine alm2map_iqu_hpx(nside, lmax, mmax, alms, Miqu)
	integer p, l, m, n, nside, lmax, mmax
	complex(DPC) alms(1,0:lmax,0:mmax)
	real(DP) Miqu(0:12*nside**2-1,3), theta, phi
	
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
	real(DP) Miqu(0:12*nside**2-1,3), theta, phi
	
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print location info (map is assumed to be in nested order)
! kind specifies peak (-1/+1 for min/max) or zero crossing (0)
subroutine location(nside, pixel, value, kind)
	integer nside, pixel, kind
	real(DP) theta, phi, value
	
	call pix2ang_nest(nside, pixel, theta, phi)
	write (*,'(3G24.16,I)') theta, phi, value, kind
end subroutine

! find extrema in the interior of the masked map
subroutine extrema(M, mask, nside, n)
	integer i, k, n, nside, nn(8)
	real(DP), dimension(0:n) :: M, mask
	
	do i = 0,n
		if (mask(i) == 0.0) cycle
		call neighbours_nest(nside, i, nn, k)
		if (any(mask(nn(1:k)) == 0.0)) cycle
		
		if (M(i) > maxval(M(nn(1:k)))) call location(nside, i, M(i), +1)
		if (M(i) < minval(M(nn(1:k)))) call location(nside, i, M(i), -1)
	end do
end subroutine

! find skeleton in the interior of the masked map
subroutine skeleton(M, mask, nside, n)
	integer i, k, n, nside, nn(8)
	real(DP), dimension(0:n) :: M, mask
	
	do i = 0,n
		if (mask(i) == 0.0) cycle
		call neighbours_nest(nside, i, nn, k)
		if (any(mask(nn(1:k)) == 0.0)) cycle
		
		if (any(M(i)*M(nn(1:k)) < 0.0)) call location(nside, i, M(i), 0)
	end do
end subroutine

! calculate cumulative Minkowski functionals
subroutine minkowski(M, mask, nmaps, n, bins)
	integer i, k, n, p, nmaps, bins, s, o
	real(DP) F(nmaps), M(0:n,nmaps), mask(0:n)
	real(SP), allocatable :: U(:)
	integer, allocatable :: idx(:), pxl(:)
	
	! temporary storage
	allocate(U(n+1), idx(n+1), pxl(n+1))
	
	! index unmasked pixels
	k = 0; do i = 0,n; if (mask(i) /= 0.0) then
		k = k+1; U(k) = M(i,1); pxl(k) = i
	end if; end do
	
	call indexx(k, U, idx)
	
	! output accumulated distribution
	F = 0.0; s = k/bins; o = (k - bins*s)/2 + 1
	
	do i = 1,k; p = pxl(idx(i))
		F = F + (/ 1.0, M(p,2:) /)
		if (mod(i-o,s) == 0) write (*,'(8G24.16)') M(p,1), F/k
	end do
	
	! clean up after ourselves
	deallocate(U, idx, pxl)
end subroutine

end
