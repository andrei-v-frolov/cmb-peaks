! $Id$
! Remap statistical estimator map into p-value from a bunch of simulations
! invoke: remap {utp|ltp|mtp} <map.fits> <out.fits> <sims-000.fits> ... <sims-XXX.fits>

program remap

! HEALPix includes
use mapio
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: mapping, fin, fout, fsim
integer :: i, m, n, nsims, nmaps = 0, nside = 0, ord = 0
real(IO), allocatable :: Map(:,:), Mout(:,:), sims(:,:,:)

integer, parameter :: LTP = 1	! lower tail probability (i.e. #sims < v)
integer, parameter :: UTP = 2	! upper tail probability (i.e. #sims > v)
integer, parameter :: MTP = 3	! multi tail probability (i.e. #sims < v for lower half, and > v for upper half)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! argument sanity check
if (nArguments() < 3) call abort("usage: remap {utp|ltp|mtp} <map.fits> <out.fits> <sims-****.fits>")
if (nArguments() < 5) call abort("Too few simulations supplied, cannot remap...")

! parse argument list
call getArgument(1, mapping)
call getArgument(2, fin )
call getArgument(3, fout)
nsims = nArguments() - 3

! read input map
call read_map(fin, Map, nside, nmaps, ord); n = nside2npix(nside)-1

! allocate dynamic arrays
allocate(Mout, mold=Map)
allocate(sims(0:n,nmaps,nsims))

! read simulated maps
do i = 1,nsims
	call getArgument(i+3, fsim)
	call read_map(fsim, Mout, nside, nmaps, ord); sims(:,:,i) = Mout
end do

! remap pixels into log-P using distribution of simulated values
select case (mapping)
	case ('utp','UTP'); call logp_map(UTP)
	case ('ltp','LTP'); call logp_map(LTP)
	case ('mtp','MTP','mUTP'); call logp_map(MTP)
	case default; call abort("Mapping not supported, use one of UTP, LTP, or mUTP...")
end select

! output log-P map
call write_map(fout, Mout, nside, ord, creator='REMAP')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! remap statistical estimator map into p-value 
subroutine logp_map(mapping)
	integer i, m, mapping
	
	do m = 1,nmaps; do i = 0,n
		Mout(i,m) = -logutp(mapping, Map(i,m), nsims, sims(i,m,:), clip=.true.)
	end do; end do
end subroutine

! calculate log of upper tail probability of v (i.e. #sims > v)
function logutp(mapping, v, n, samples, bootstrap, clip)
	integer mapping, n; real(IO) logutp, v, samples(n)
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
	select case (mapping)
		case (LTP); forall (i=1:m) F(i) = log10(real(i)/real(m))
		case (UTP); forall (i=1:m) F(i) = log10(real(m+1-i)/real(m))
		case (MTP); forall (i=1:m) F(i) = log10(min(i,m+1-i)/real(m))
		case default; call abort("Mapping not supported...")
	end select
	
	! linear interpolation of log-UTP
	i = count(X < v); if (i < 1) i = 1; if (i > m-1) i = m-1
	logutp = F(i) + (F(i+1)-F(i)) * (v-X(i))/(X(i+1)-X(i))
	
	! clip if requested
	if (present(clip) .and. clip .and. v < X(1)) logutp = F(1)
	if (present(clip) .and. clip .and. v > X(m)) logutp = F(m)
end function

end
