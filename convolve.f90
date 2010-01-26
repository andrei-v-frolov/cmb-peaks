! $Id$
! brute-force convolution of a (masked) map with a circular kernel
! invoke: convolve <beam fwhm in acrmin> <map.fits> <out.fits[:downgrade]> [mask.fits[:tolerance]]

program convolve

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: dg = 3; real :: eps = 1.0
character(len=80) :: header(64), beam, fin, fout, fmask
integer nmaps, nside, npix, n, mside, mpix, m, ntot, ord
real(SP), allocatable :: Min(:,:), Mout(:,:), Mask(:,:)

integer i, j, k, l
integer, allocatable :: idx(:)
real(DP) :: vi(3), vj(3), t, p, aw, ww
real(DP), allocatable :: w(:), s(:), q(:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! beam FWHM
call getArgument(1, beam);
read (beam,*) t; t = t/(180.0*60.0) * pi
ww = 2.0*sin(t/2.0)/fwhm(); aw = 2.0*asin(ww/2.0)

! input map
call getArgument(2, fin )

! output map and optional downgrade exponent
call getArgument(3, fout); i = index(fout, ":", .true.)
if (i > 0) then; read (fout(i+1:),*) dg; fout(i:) = ""; end if

! import input map
ntot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)

mside = nside/2**dg
mpix = nside2npix(mside); m = mpix-1
npix = nside2npix(nside); n = npix-1

allocate(idx(0:n), w(nmaps), s(nmaps), q(nmaps))
allocate(Min(0:n,nmaps), Mout(0:m,nmaps), Mask(0:n,nmaps))

call input_map(fin, Min, npix, nmaps)
if (ord /= 1) call convert_nest2ring(nside, Min)

! import mask if specified
if (nArguments() < 4) then; Mask = 1.0; else
	call getArgument(4, fmask); i = index(fmask, ":", .true.)
	if (i > 0) then; read (fmask(i+1:),*) eps; fmask(i:) = ""; end if
	
	ntot = getsize_fits(fmask, ordering=ord)
	
	call input_map(fmask, Mask, npix, nmaps)
	if (ord /= 1) call convert_nest2ring(nside, Mask)
end if

! convolution loop
do j = 0,m
	s = 0.0 ! value  accumulator
	q = 0.0 ! weight accumulator
	p = 0.0 ! kernel accumulator
	
	call pix2vec_ring(mside, j, vj)
	call query_disc(nside, vj, aw, idx, k)
	
	do l = 0,k-1; i = idx(l)
		call pix2vec_ring(nside, i, vi)
		call angdist(vi, vj, t)
		
		t = 2.0*sin(t/2.0)
		t = kernel(t/ww)
		w = t*Mask(i,:)
		
		s = s + w*Min(i,:)
		q = q + w
		p = p + t
	end do
	
	! mask away non-conforming pixels
	where (abs((q-p)/p) > eps/2.0) q = 0.0
	
	Mout(j,:) = s/q
end do

call write_minimal_header(header, 'MAP', nside=mside, order=1)

call output_map(Mout, header, '!'//fout)

contains

! convolution kernel
function kernel(t)
	real t, x, kernel
	
	! compact support
	if (t >= 1.0) then
		kernel = 0.0
		return
	end if
	
	! symmetric kernel
	x = t*t
	
	! polynomial filters (generalized Savitzky-Golay family, interval orthogonality)
	!kernel = (1.0/2.0)																	! SG00
	!kernel = (3.0/4.0)*(1.0-x)																! SG01
	!kernel = (15.0/16.0)*(1.0-x)**2															! SG02
	!kernel = (225.0/128.0+(-525.0/64.0+945.0/128.0*x)*x)													! SG40
	!kernel = (525.0/256.0+(-1575.0/128.0+3465.0/256.0*x)*x)*(1.0-x)											! SG41
	!kernel = (4725.0/2048.0+(-17325.0/1024.0+45045.0/2048.0*x)*x)*(1.0-x)**2										! SG42
	!kernel = (99225.0/32768.0+(-363825.0/8192.0+(2837835.0/16384.0+(-2027025.0/8192.0+3828825.0/32768.0*x)*x)*x)*x)					! SG80
	!kernel = (218295.0/65536.0+(-945945.0/16384.0+(8513505.0/32768.0+(-6891885.0/16384.0+14549535.0/65536.0*x)*x)*x)*x)*(1.0-x)				! SG81
	!kernel = (945945.0/262144.0+(-4729725.0/65536.0+(48243195.0/131072.0+(-43648605.0/65536.0+101846745.0/262144.0*x)*x)*x)*x)*(1.0-x)**2			! SG82
	!kernel = (34459425.0/8388608.0+(-218243025.0/2097152.0+(2749862115.0/4194304.0+(-3011753745.0/2097152.0+8365982625.0/8388608.0*x)*x)*x)*x)*(1.0-x)**4	! SG84
	
	! polynomial filters (generalized Savitzky-Golay family, spherical orthogonality)
	!kernel = 2.0																		! SSG00
	!kernel = 4.0*(1.0-x)																	! SSG01
	 kernel = 6.0*(1.0-x)**2																! SSG02
	!kernel = (18.0+(-72.0+60.0*x)*x)															! SSG40
	!kernel = (24.0+(-120.0+120.0*x)*x)*(1.0-x)														! SSG41
	!kernel = (30.0+(-180.0+210.0*x)*x)*(1.0-x)**2														! SSG42
	!kernel = (50.0+(-600.0+(2100.0+(-2800.0+1260.0*x)*x)*x)*x)												! SSG80
	!kernel = (60.0+(-840.0+(3360.0+(-5040.0+2520.0*x)*x)*x)*x)*(1.0-x)											! SSG81
	!kernel = (70.0+(-1120.0+(5040.0+(-8400.0+4620.0*x)*x)*x)*x)*(1.0-x)**2											! SSG82
	!kernel = (90.0+(-1800.0+(9900.0+(-19800.0+12870.0*x)*x)*x)*x)*(1.0-x)**4										! SSG84
end function kernel

! full width half maximum of the kernel
function fwhm()
	integer i; real x, a, b, fa, fb, f0, fwhm
	
	a = 0.0;     fa = kernel(a)
	b = 1.0/3.0; fb = kernel(b)
	
	f0 = fa/2.0
	
	do i = 1,32
		if (a == b) exit
		x = b - (b-a)/(fb-fa) * (fb-f0)
		a = b; fa = fb; b = x; fb = kernel(x)
	end do
	
	fwhm = 2.0*x
end function fwhm
end
