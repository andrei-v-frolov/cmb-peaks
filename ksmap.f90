! $Id$
! Kolmogorov-Smirnov deviation of peaks in a circular aperture from full sky
! invoke: cat peaks.dat | ksmap <out.fits> <disk diameter, deg>

program ksmap

! HEALPix includes
use mapio
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: pks = 12 * 512**2 ! max number of peaks we can handle

real diam, qcut
real(DP) peak(4,pks)
integer i, n, m, status

character(len=80) :: fout, fwhm
integer :: nside = 64, ord = NEST
real(IO), allocatable :: Map(:,:)
real(DP), allocatable :: V(:,:)
integer, allocatable :: idx(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse argument list
if (nArguments() /= 2) call abort("usage: cat peaks.dat | ksmap <out.fits> [disk diameter, deg]")

call getArgument(1, fout)
call getArgument(2, fwhm)

read (fwhm,*) diam; qcut = cos(diam*pi/360.0)

! read peak data from stdin
do i = 1,pks
	read (*,*,iostat=status) peak(:,i)
	if (status < 0) exit
end do; n = i-1

! check if we overflowed by any chance
if (i > pks) call abort("ksmap: Peak count limit exceeded, recompile with larger value for pks...")

! allocate storage
allocate(Map(0:12*nside**2-1,2), V(3,n), idx(n))

! sort the distribution
call indexx(n, peak(3,1:n), idx)

! peak positions, sorted in increasing temperature order
do i = 1,n
	call ang2vec(peak(1,idx(i)), peak(2,idx(i)), V(:,i))
end do

! scan disk aperture around the sky
do i = 1,12*nside**2-1
	call select_peaks(i, qcut, idx, m, n)
	Map(i,1) = ksdev(idx, m, n); Map(i,2) = m
end do

! output Kolmogorov-Smirnov deviation map
call write_map(fout, Map, nside, ord, creator='KSMAP')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! select peaks within a disk centered around pixel
subroutine select_peaks(pxl, qcut, idx, m, n)
	real qcut, O(3); integer i, k, m, n, idx(n), pxl
	
	call pix2vec_nest(nside, pxl, O)
	
	k = 1; do i = 1,n
		if (sum(V(:,i)*O) < qcut) cycle
		idx(k) = i; k = k+1
	end do; m = k-1
end subroutine

! calculate Kolmogorov-Smirnov deviation of subsample
function ksdev(idx, m, n)
	integer m, n, idx(n); real ksdev
	
	ksdev = sqrt(real(n*m)/real(n+m)) * maxval(abs( ((/1:m/)-1.0)/(m-1.0) - (idx(1:m)-1.0)/(n-1.0) ))
end function ksdev

end
