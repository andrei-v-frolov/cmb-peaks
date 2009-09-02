! $Id$
! convolve  -  brute-force convolution of a (masked) map

program convolve

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: header(64), fin, fout, fmask
integer nmaps, nside, npix, n, mside, mpix, m, ntot, ord
real(SP), allocatable :: Min(:,:), Mout(:,:), Mask(:,:)

integer i, j, k, l
real(DP) :: vi(3), vj(3), t
integer, allocatable :: idx(:)
real(DP), allocatable :: w(:), s(:), q(:)

integer, parameter :: dg = 4
real, parameter :: ww = 5.0/180.0 * pi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call getArgument(1, fin)
call getArgument(2, fout)

ntot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)

mside = nside/2**dg
mpix = nside2npix(mside); m = mpix-1
npix = nside2npix(nside); n = npix-1

allocate(idx(0:n), w(nmaps), s(nmaps), q(nmaps))
allocate(Min(0:n,nmaps), Mout(0:m,nmaps), Mask(0:n,nmaps))

call input_map(fin, Min, npix, nmaps)

! import mask if specified
if (nArguments() > 2) then
	call getArgument(3, fmask)
	call input_map(fmask, Mask, npix, nmaps)
else
	Mask = 1.0
end if

! convolution loop
do j = 0,m
	s = 0.0 ! value  accumulator
	q = 0.0 ! weight accumulator
	
	call pix2vec(mside, j, vj)
	call query_disc(nside, vj, ww, idx, k, ord-1)
	
	do l = 0,k-1; i = idx(l)
		call pix2vec(nside, i, vi)
		call angdist(vi, vj, t)
		
		w = kernel(t/ww)*Mask(i,:)
		
		s = s + w*Min(i,:)
		q = q + w
	end do
	
	Mout(j,:) = s/q
	
	write (*,*) j,k,s/q
end do

call write_minimal_header(header, 'MAP', nside=mside, order=ord)

call output_map(Mout, header, '!'//fout)

contains

! convolution kernel
function kernel(x)
	real x, kernel; kernel = 0.0
	
	if (x < 1.0) kernel = (1.0 - x**2)**2
end function kernel

! renders cartesian vector coordinates of the nominal pixel center
subroutine pix2vec(n, p, v)
	integer n, p; real(DP) :: v(3); v = 0.0
	
	if (ord == 1) call pix2vec_ring(n, p, v)
	if (ord == 2) call pix2vec_nest(n, p, v)
end subroutine pix2vec

end
