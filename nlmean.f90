! HEALPix non-local mean filter - noise cleaning for non-Gaussian maps
! usage: nlmean FWHM amount map.fits output.fits [residual.fits]

program nlmeanf

! HEALPix includes
use mapio
use extension
use almtools

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real, parameter :: eps = 1.0e-7

character(len=8000) :: arg, fin, fout, fres
integer :: nmaps = 0, nside = 0, ord = 0
real(IO), dimension(:,:), allocatable :: Minp, Mout
real(DP), dimension(:,:), allocatable :: M, F, S
integer, dimension(:,:), allocatable :: idx, rank

real fwhm, amount
real(DP) sigma(3)
integer i, n, bmax, lmax, status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nArguments() < 4 .or. nArguments() > 5) call abort("sage: nlmean FWHM amount map.fits output.fits [residual.fits]")

call getArgument(1, arg); read(arg, *, iostat=status) fwhm; if (status /= 0) call abort("error parsing FWHM specification: " // trim(arg))
call getArgument(2, arg); read(arg, *, iostat=status) amount; if (status /= 0) call abort("error parsing filtering amount: " // trim(arg))
call getArgument(3, fin); if (verbose) write (*,*) "Reading " // trim(fin); call read_map(fin, Minp, nside, nmaps, ord)
call getArgument(4, fout); if (nArguments() == 5) call getArgument(5, fres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bmax = int(180*240*sqrt(-log(2.0)*log(eps))/(pi*fwhm) + 0.5)
n = nside2npix(nside)-1; lmax = min(3*nside-1, bmax, 1000)

! allocate workspace and output storage
allocate(Mout, mold=Minp); Mout = 0.0
allocate(M(0:n,nmaps), S(0:n,nmaps), F(0:n,3), idx(n+1,3), rank(n+1,3))

! bring input map into RING ordering (if required)
if (ord == NEST) call convert_nest2ring(nside, Minp); M = Minp

! extract features and build rank tables
if (verbose) write (*,*) "Extracting features from the map..."
call features(nside, lmax, fwhm, M(:,1), F)

if (verbose) write (*,*) "Ranking feature maps and estimating variance..."
do i = 1,3
	call indexx(n+1, F(:,i), idx(:,i))
	call rankx(n+1, idx(:,i), rank(:,i))
end do

! estimate feature variance using L2-moment...
sigma = sum(F*(2*rank-n-1), 1)/(n-1.0)**2 * sqrt(pi)
write (*,*) sigma

! ...
if (verbose) write (*,*) "Applying brute-force non-local mean filter..."
call nlbrute(nside, nmaps, 1.0/(amount*sigma)**2, F, M, S)

! write output map(s)
if (verbose) write (*,*) "Saving filtered map to " // trim(fout)
Mout = S; if (ord == NEST) call convert_ring2nest(nside, Mout)
call write_map(fout, Mout, nside, ord, pol=0, vec=0, creator='NLMEAN')

contains

! synthesize smoothed feature map based on Minkowski functionals
subroutine features(nside, lmax, fwhm, map, F)
	integer nside, lmax, l; real fwhm
	real(DP) map(0:12*nside**2), F(0:12*nside**2,3), bls(0:lmax,1)
	complex(DPC), allocatable :: alms(:,:,:)
	
	allocate(alms(1,0:lmax,0:lmax))
	
	call map2alm(nside, lmax, lmax, map, alms)
	call gaussbeam(fwhm, lmax, bls)
	
	forall (l=0:lmax) alms(1,l,0:l) = bls(l,1)*alms(1,l,0:l)
	call alm2map_covariant(nside, lmax, lmax, alms, mink=F)
	
	deallocate(alms)
end subroutine

! Gaussian weight kernel evaluating similarity of feature indicators
pure function weight(w, v)
	real(DP) weight, w(3), v(3); intent(in) w, v
	
	weight = exp(-sum(w*v*v)/2.0)
end function

! ...
pure function nlmean(nside, nmaps, w, v, F, map)
	integer nside, nmaps, i, n
	real(DP) w(3), v(3), nlmean(nmaps), s(0:nmaps)
	real(DP), dimension(0:12*nside**2, 3) :: F
	real(DP), dimension(0:12*nside**2, nmaps) :: map
	intent(in) nside, nmaps, w, v, F, map
	
	n = 12*nside**2-1; s = 0.0
	
	do i = 0,n
		s = s + weight(w, v-F(i,:)) * [1.0, map(i,:)]
	end do
	
	nlmean = s(1:nmaps)/s(0)
end function

! brute-force non-local means filter
subroutine nlbrute(nside, nmaps, w, F, map, out)
	integer nside, nmaps, i, n; real(DP) w(3)
	real(DP), dimension(0:12*nside**2, 3) :: F
	real(DP), dimension(0:12*nside**2, nmaps) :: map, out
	intent(in) nside, nmaps, w, F, map; intent(out) out
	
	n = nside2npix(nside)-1
	
	!$OMP PARALLEL DO
	do i = 0,n; out(i,:) = nlmean(nside, nmaps, w, F(i,:), F, map); end do
	!$OMP END PARALLEL DO
end subroutine

end