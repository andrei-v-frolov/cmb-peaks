! $Id$
! Remap statistical estimator map into p-value from a bunch of simulations
! invoke: remap <map.fits> <out.fits> <sims-000.fits> ... <sims-XXX.fits>

program remap

! HEALPix includes
use mapio
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: fin, fout, fsim
integer :: i, m, n, nsims, nmaps = 0, nside = 0, ord = 0
real(IO), allocatable :: Min(:,:), Mout(:,:), sims(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse argument list
call getArgument(1, fin )
call getArgument(2, fout)
nsims = nArguments() - 2

if (nsims < 1) call abort("Too few simulations supplied, cannot remap...")

! read input map
call read_map(fin, Min, nside, nmaps, ord); n = nside2npix(nside)-1

! allocate dynamic arrays
allocate(Mout, mold=Min)
allocate(sims(0:n,nmaps,nsims))

! read simulated maps
do i = 1,nsims
	call getArgument(i+2, fsim)
	call read_map(fsim, Mout, nside, nmaps, ord); sims(:,:,i) = Mout
end do

! remap pixels according 
do m = 1,nmaps; do i = 0,n
	Mout(i,m) = logutp(Min(i,m), nsims, sims(i,m,:), clip=.true.)
end do; end do


! output L-weight masks (to a single FITS container)
call write_map(fout, Mout, nside, ord, creator='REMAP')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function logutp(v, n, samples, bootstrap, clip)
	integer n; real(IO) logutp, v, samples(n)
	logical, optional, value :: bootstrap, clip
	
	! local storage
	real(DP) X(n), F(n)
	integer i, m, idx(n)
	logical valid(n)
	
	! pack valid samples
	valid = .not. isnan(samples)
	X = pack(samples, valid)
	m = count(valid)
	
	! optionally bootstrap
	if (present(bootstrap) .and. bootstrap) then
		call random_number(F)
		idx = floor(F*m) + 1
		X(1:m) = X(idx(1:m))
	end if
	
	! sort the distribution
	call indexx(m, X, idx)
	X(1:m) = X(idx(1:m))
	
	! form a log-UTP for samples
	forall (i=1:m) F(i) = log10((m+1-i)/real(m))
	
	! linear interpolation of log-UTP
	i = count(X < v); if (i < 1) i = 1; if (i > m-1) i = m-1
	logutp = F(i) + (F(i+1)-F(i)) * (v-X(i))/(X(i+1)-X(i))
	
	! clip if requested
	if (present(clip) .and. clip .and. v < X(1)) logutp = F(1)
	if (present(clip) .and. clip .and. v > X(m)) logutp = F(m)
end function

end
