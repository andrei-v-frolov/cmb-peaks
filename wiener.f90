! $Id$
! Generate Wiener-optimal filter of a circular kernel subject to noise with known Cl's
! invoke: wiener <{kernel:fwhm|kernel-bls.fits}> [noise-cls.{ascii|fits}] [filter-bls.fits]

program wiener

! HEALPix includes
use extension
use head_fits
use healpix_types
use fitstools

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: lmax = 4000, ncl = 1
character(len=80) :: header(64), kern, fin, fout
real(DP), dimension(0:lmax,ncl) :: Bl, Cl

integer i, l; real t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! filter kernel
call getArgument(1, kern)
i = index(kern, ":", .true.)

if (i > 0) then
	read (kern(i+1:),*) t; kern(i:) = ""
	call beam(Bl, lmax, kern, t * pi/(180.0*60.0))
else
	call data2cl(kern, Bl, lmax, ncl, header)
end if

! noise spectrum
if (nArguments() < 2) then; Cl = 1.0; else
	call getArgument(2, fin)
	call data2cl(fin, Cl, lmax, ncl, header)
	
	! fill in missing multipoles
	call intcls(Cl(0:8,1), 8)
end if

! Wiener filter
where (Cl /= 0.0) Bl = Bl/sqrt(Cl)

! output filter
if (nArguments() < 3) then
	do l = 0,lmax; write (*,'(i,8G24.16)') l, Bl(l,:); end do
else
	call getArgument(3, fout)
	
	call write_minimal_header(header, 'CL', nlmax=lmax)
	call write_asctab(Bl, lmax, ncl, header, 64, '!'//fout)
end if


contains

! calculate (scaled) Legendre polynomials up to order lmax at x
subroutine legendre(P, lmax, x, s)
	integer l, lmax; real(DP) P(0:lmax), x, s
	
	P(0) = s
	P(1) = s*x
	
	do l = 2,lmax
		P(l) = ((2*l-1)*x*P(l-1) - (l-1)*P(l-2))/l
	end do
end subroutine

! calculate beam bl's for a kernel of angular FWHM alpha
subroutine beam(Bl, lmax, k, alpha)
	integer, parameter :: nmin = 8, nmax = 13
	integer lmax; character(*) k; real x, w, alpha
	real(DP) Bl(0:lmax), P(0:lmax), S(0:nmax,0:lmax), h(0:nmax)
	integer i, ii, l, n
	
	! zero width means delta beam
	if (alpha == 0.0) then
		Bl = 1.0
		return
	end if
	
	! assuming FIR kernel on unit disk
	w = 2.0*sin(alpha/2.0)/fwhm(k); h(0) = 1.0
	call legendre(S(0,:), lmax, 1.0-w**2/2.0, kernel(k,1.0)/2.0)
	
	! Romberg integration loop
	do n = 1,nmax; ii=2**(n-1)
		! initialize from coarse grid
		S(n,:) = S(n-1,:)/2.0; h(n) = h(n-1)/4.0
		
		! accumulate on refined points
		do i = 1,ii; x = (i-0.5)/ii
			call legendre(P, lmax, 1.0-(w*x)**2/2.0, kernel(k,x)*x/ii/2.0); S(n,:) = S(n,:) + P
		end do
	end do
	
	! interpolate to zero stepsize
	do l=0,lmax
		call polint(h(nmin:nmax), S(nmin:nmax,l), nmax-nmin+1, 0.0, Bl(l), x)
	end do
end subroutine

! convolution kernel
function kernel(k, t)
	real t, x, kernel; character(*) k
	
	! compact support
	if (t > 1.0) then
		kernel = 0.0
		return
	end if
	
	! symmetric kernel
	x = t*t
	
	select case (k)
		! polynomial filters (generalized Savitzky-Golay family, spherical orthogonality)
		case('SSG00');	kernel = 2.0
		case('SSG01');	kernel = 4.0-4.0*x
		case('SSG02');	kernel = 6.0*(1.0-x)**2
		case('SSG04');	kernel = 10.0*(1.0-x)**4
		case('SSG08');	kernel = 18.0*(1.0-x)**8
		case('SSG20');	kernel = 8.0-12.0*x
		case('SSG21');	kernel = (12.0-24.0*x)*(1.0-x)
		case('SSG22');	kernel = (16.0-40.0*x)*(1.0-x)**2
		case('SSG24');	kernel = (24.0-84.0*x)*(1.0-x)**4
		case('SSG28');	kernel = (40.0-220.0*x)*(1.0-x)**8
		case('SSG40');	kernel = 18.0+(-72.0+60.0*x)*x
		case('SSG41');	kernel = (24.0+(-120.0+120.0*x)*x)*(1.0-x)
		case('SSG42');	kernel = (30.0+(-180.0+210.0*x)*x)*(1.0-x)**2
		case('SSG44');	kernel = (42.0+(-336.0+504.0*x)*x)*(1.0-x)**4
		case('SSG48');	kernel = (66.0+(-792.0+1716.0*x)*x)*(1.0-x)**8
		case('SSG80');	kernel = 50.0+(-600.0+(2100.0+(-2800.0+1260.0*x)*x)*x)*x
		case('SSG81');	kernel = (60.0+(-840.0+(3360.0+(-5040.0+2520.0*x)*x)*x)*x)*(1.0-x)
		case('SSG82');	kernel = (70.0+(-1120.0+(5040.0+(-8400.0+4620.0*x)*x)*x)*x)*(1.0-x)**2
		case('SSG84');	kernel = (90.0+(-1800.0+(9900.0+(-19800.0+12870.0*x)*x)*x)*x)*(1.0-x)**4
		case('SSG88');	kernel = (130.0+(-3640.0+(27300.0+(-72800.0+61880.0*x)*x)*x)*x)*(1.0-x)**8
		case('SSGX0');	kernel = 162.0+(-6480.0+(83160.0+(-498960.0+(1621620.0+(-3027024.0&
						+(3243240.0+(-1853280.0+437580.0*x)*x)*x)*x)*x)*x)*x)*x
		case('SSGX1');	kernel = (180.0+(-7920.0+(110880.0+(-720720.0+(2522520.0+(-5045040.0&
						+(5765760.0+(-3500640.0+875160.0*x)*x)*x)*x)*x)*x)*x)*x)*(1.0-x)
		case('SSGX2');	kernel = (198.0+(-9504.0+(144144.0+(-1009008.0+(3783780.0+(-8072064.0&
						+(9801792.0+(-6301152.0+1662804.0*x)*x)*x)*x)*x)*x)*x)*x)*(1.0-x)**2
		case('SSGX4');	kernel = (234.0+(-13104.0+(229320.0+(-1834560.0+(7796880.0+(-18712512.0&
						+(25395552.0+(-18139680.0+5290740.0*x)*x)*x)*x)*x)*x)*x)*x)*(1.0-x)**4
		case('SSGX8');	kernel = (306.0+(-22032.0+(488376.0+(-4883760.0+(25639740.0+(-75209904.0&
						+(123559128.0+(-105907824.0+36773550.0*x)*x)*x)*x)*x)*x)*x)*x)*(1.0-x)**8
		
		! (truncated) Gaussian kernel
		case('GAUSS');		kernel = 50.0*exp(-25.0*x)
	end select
end function kernel

! full width half maximum of the kernel
function fwhm(k)
	integer i; real x, a, b, fa, fb, f0, fwhm; character(*) k
	
	! kernel maximum is at origin
	a = 0.0; fa = kernel(k,a); f0 = fa/2.0
	
	! hunt for the half maximum
	do i = 1,32
		b = i/32.0; fb = kernel(k,b)
		if (fb < f0) exit
		a = b; fa = fb
	end do
	
	! refine using Newton's method
	do i = 1,32
		if (a == b) exit
		x = b - (b-a)/(fb-fa) * (fb-f0)
		a = b; fa = fb; b = x; fb = kernel(k,x)
	end do
	
	fwhm = 2.0*x
end function fwhm

! read (normalized) Cl's from a data file
subroutine data2cl(fin, Cl, lmax, ncl, header)
	integer l, lmax, ncl, status
        character(*) fin, header(:)
	real(DP) p(ncl), Cl(0:lmax,ncl)
	
	! if it looks like a FITS file, read it in unchanged
	if (index(fin, ".fits") > 0) then
		call fits2cl(fin, Cl, lmax, ncl, header); return
	end if
	
	! nope, it's a plain text file (with normalized Cl's!)
        Cl = 0.0; open(11, file=fin, action="read")
        
        do
		read (11,*,iostat=status) l, p
		if (status < 0) exit
		if (l > lmax) cycle
		Cl(l,:) = 2.0*pi*p/(l*(l+1.0))
	end do
        
        close(11)
end subroutine

! interpolate missing lower multipoles
subroutine intcls(Cl, lmax)
	integer i, l, lmax; real dy
	real(DP), dimension(0:lmax) :: Cl, x, y, w
	
	! interpolation is done on a scaled Cl's
	x = (/0:lmax/); w = x*x+0.7071067812; y = w*Cl
	
	! find what multipoles are missing
	do l = 0,lmax/2; if (Cl(l) > 0.0) exit; end do
	
	! interpolate missing multipoles
	do i = 0,l-1
		call polint(x(l:lmax), y(l:lmax), lmax-l+1, x(i), y(i), dy); Cl(i)=y(i)/w(i)
	end do
end subroutine

end
