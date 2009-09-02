! $Id$
! brute-force convolution of a (masked) map with a circular kernel
! invoke: convolve ksize(acrmin) map.fits out.fits [mask.fits]

program convolve

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: header(64), win, fin, fout, fmask
integer nmaps, nside, npix, n, mside, mpix, m, ntot, ord
real(SP), allocatable :: Min(:,:), Mout(:,:), Mask(:,:)

integer i, j, k, l
integer, allocatable :: idx(:)
real(DP) :: vi(3), vj(3), t, p, ww
real(DP), allocatable :: w(:), s(:), q(:)

integer, parameter :: dg = 4		! downgrade output
real, parameter :: eps = 0.10		! masking tolerance


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call getArgument(1, win)
call getArgument(2, fin)
call getArgument(3, fout)

read (win,*) t; ww = t/(180.0*60.0) * pi

ntot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)

mside = nside/2**dg
mpix = nside2npix(mside); m = mpix-1
npix = nside2npix(nside); n = npix-1

allocate(idx(0:n), w(nmaps), s(nmaps), q(nmaps))
allocate(Min(0:n,nmaps), Mout(0:m,nmaps), Mask(0:n,nmaps))

call input_map(fin, Min, npix, nmaps)
if (ord /= 1) call convert_nest2ring(nside, Min)

! import mask if specified
if (nArguments() > 3) then
	call getArgument(4, fmask)
	
	ntot = getsize_fits(fmask, ordering=ord)
	
	call input_map(fmask, Mask, npix, nmaps)
	if (ord /= 1) call convert_nest2ring(nside, Mask)
else
	Mask = 1.0
end if

! convolution loop
do j = 0,m
	s = 0.0 ! value  accumulator
	q = 0.0 ! weight accumulator
	p = 0.0 ! kernel accumulator
	
	call pix2vec_ring(mside, j, vj)
	call query_disc(nside, vj, ww, idx, k)
	
	do l = 0,k-1; i = idx(l)
		call pix2vec_ring(nside, i, vi)
		call angdist(vi, vj, t)
		
		t = kernel(t/ww)
		w = t*Mask(i,:)
		
		s = s + w*Min(i,:)
		q = q + w
		p = p + t
	end do
	
	! mask away non-conforming pixels
	where (abs((q-p)/p) > eps) q = 0.0
	
	Mout(j,:) = s/q
end do

call write_minimal_header(header, 'MAP', nside=mside, order=1)

call output_map(Mout, header, '!'//fout)

contains

! convolution kernel
function kernel(x)
	real x, kernel; kernel = 0.0
	
	if (x < 1.0) kernel = (1.0 - x**2)**2
end function kernel

end
