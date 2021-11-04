! HEALPix non-local mean filter - noise cleaning for non-Gaussian maps
! usage: nlmean FWHM amount map.fits output.fits [residual.fits]

program nlmean

! HEALPix includes
use mapio
use extension
use almtools

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=8000) :: arg, fin, fout, fres
integer :: nmaps = 0, nside = 0, ord = 0
real(IO), dimension(:,:), allocatable :: Min, Mout
real(DP), dimension(:,:), allocatable :: M, F

real fwhm, amount
integer n, lmax, status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nArguments() < 4 .or. nArguments() > 5) call abort("sage: nlmean FWHM amount map.fits output.fits [residual.fits]")

call getArgument(1, arg); read(arg, *, iostat=status) fwhm; if (status /= 0) call abort("error parsing FWHM specification: " // trim(arg))
call getArgument(2, arg); read(arg, *, iostat=status) amount; if (status /= 0) call abort("error parsing filtering amount: " // trim(arg))
call getArgument(3, fin); call read_map(fin, Min, nside, nmaps, ord)
call getArgument(4, fout); if (nArguments() == 5) call getArgument(5, fres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n = nside2npix(nside)-1; lmax = 3*nside-1; lmax=1000

! allocate workspace and output storage
allocate(Mout, mold=Min); Mout = 0.0
allocate(M(0:n,nmaps), F(0:n,3))

if (ord == NEST) call convert_nest2ring(nside, Min); M = Min

call features(nside, lmax, fwhm, M(:,1), F)

Mout = F(:,1:nmaps); if (ord == NEST) call convert_ring2nest(nside, Mout)

! write output map(s)
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

end