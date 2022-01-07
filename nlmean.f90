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

real(DP) fwhm, amount, sigma, w(3)
integer i, n, bmax, lmax, status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nArguments() < 4 .or. nArguments() > 5) call abort("Usage: nlmean FWHM amount map.fits output.fits [residual.fits]")

call getArgument(1, arg); read(arg, *, iostat=status) fwhm; if (status /= 0) call abort("error parsing FWHM specification: " // trim(arg))
call getArgument(2, arg); read(arg, *, iostat=status) amount; if (status /= 0) call abort("error parsing filtering amount: " // trim(arg))
call getArgument(3, fin); if (verbose) write (*,*) "Reading " // trim(fin); call read_map(fin, Minp, nside, nmaps, ord)
call getArgument(4, fout); if (nArguments() == 5) call getArgument(5, fres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bmax = int(180*240*sqrt(-log(2.0)*log(eps))/(pi*fwhm) + 0.5)
n = nside2npix(nside)-1; lmax = min(3*nside-1, bmax, 1000)

! allocate workspace and output storage
allocate(Mout, mold=Minp); Mout = 0.0
allocate(M(0:n,nmaps), S(0:n,nmaps), F(0:n,3))

! bring input map into RING ordering (if required)
if (ord == NEST) call convert_nest2ring(nside, Minp); M = Minp

! extract features and estimate feature covariance
if (verbose) write (*,*) "Extracting features from the map..."
call features(nside, lmax, fwhm, M(:,1), F, w)
sigma = amount*norm2(M(:,1)-F(:,1))/(n+1)

! bring all maps into NEST ordering
call convert_ring2nest(nside, M)
call convert_ring2nest(nside, F)

! ...
if (verbose) write (*,*) "Applying brute-force non-local mean filter..."
call nldisc(nside, nmaps, 0.5, w/sigma**2, F, M, S)

! write output map(s)
if (verbose) write (*,*) "Saving filtered map to " // trim(fout)
Mout = S; if (ord == RING) call convert_nest2ring(nside, Mout)
call write_map(fout, Mout, nside, ord, pol=0, vec=0, creator='NLMEAN')

! write residual map if requested
if (nArguments() == 5) then
	if (verbose) write (*,*) "Saving residual map to " // trim(fres)
	Mout = M - S; if (ord == RING) call convert_nest2ring(nside, Mout)
	call write_map(fres, Mout, nside, ord, pol=0, vec=0, creator='NLMEAN')
end if

contains

! synthesize smoothed feature map based on Minkowski functionals
subroutine features(nside, lmax, fwhm, map, F, w)
	integer nside, lmax, l, n; real fwhm
	real(DP) map(0:12*nside**2-1), F(0:12*nside**2-1,3), bls(0:lmax,1)
	real(DP) w(3), delta, rho; optional w
	complex(DPC), allocatable :: alms(:,:,:)
	real(DP), allocatable :: v(:)
	
	n = nside2npix(nside)-1
	allocate(alms(1,0:lmax,0:lmax), v(0:n))
	
	call map2alm(nside, lmax, lmax, map, alms)
	call gaussbeam(fwhm, lmax, bls)
	
	forall (l=0:lmax) alms(1,l,0:l) = bls(l,1)*alms(1,l,0:l)
	call alm2map_covariant(nside, lmax, lmax, alms, feat=F, svar=v)
	
	if (present(w)) then
		delta = fwhm*pi/180/60/sqrt(8.0*log(2.0)); rho = sum(v)/(n+1)
		w = [1.0, 2.0*delta**2, 12.0*delta**4/(3.0 + (6.0*rho-1.0)*delta**2)]
	end if
	
	deallocate(alms, v)
end subroutine

! non-local means filter with running aperture (NESTED ordering)
subroutine nldisc(nside, nmaps, radius, w, F, map, out)
	integer nside, nmaps, i, j, k, l, n
	real(DP) radius, z, v(3), w(3), s(0:nmaps)
	real(DP), dimension(0:12*nside**2-1, 3) :: F
	real(DP), dimension(0:12*nside**2-1, nmaps) :: map, out
	intent(in) nside, nmaps, radius, w, F, map; intent(out) out
	
	integer, allocatable :: idx(:)
	
	n = nside2npix(nside)-1
	l = nside2npix(nside) * sin(radius/2.0)**2
	
	allocate(idx(0:3*l/2))
	
	!$OMP PARALLEL DO PRIVATE(z,v,s,k,idx)
	do i = 0,n
		call pix2vec_nest(nside, i, v)
		call query_disc(nside, v, radius, idx, k, nest=1)
		
		s = 0.0
		
		do j = 0,k-1
			v = F(i,:)-F(idx(j),:); z = sum(w*v*v)/2.0; if (z > 5.0) cycle
			s = s + exp(-z) * [1.0, map(idx(j),:)]
		end do
		
		out(i,:) = s(1:nmaps)/s(0)
	end do
	!$OMP END PARALLEL DO
	
	deallocate(idx)
	
end subroutine

end