! $Id$
! Synthesize a rank-ordered weigth masks for calculating L-moments
! invoke: lmask <map.fits[:channel]> <lmask.fits[:moments]> [mask.fits] [lmap-base]

program lmask

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: ch = 1, nmoms = 4
integer i, nmaps, nside, npix, ntot, nuse, ord
character(len=80) :: header(64), fin, fout, fmask
real(SP), allocatable :: Min(:,:), Mout(:,:), Mask(:,:)
integer, allocatable :: indx(:), rank(:)
real(DP), allocatable :: P(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! input map and optional signal channel selection
call getArgument(1, fin ); i = index(fin,  ":", .true.)
if (i > 0) then; read (fin(i+1:),*) ch; fin(i:) = ""; end if

! output map and optional number of moments selection
call getArgument(2, fout); i = index(fout, ":", .true.)
if (i > 0) then; read (fout(i+1:),*) nmoms; fout(i:) = ""; end if

! allocate dynamic arrays to store maps and ranks
ntot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord); npix = nside2npix(nside)
allocate(Min(npix,nmaps), Mask(npix,nmaps), Mout(npix,nmoms), P(nmoms,nmoms), indx(npix), rank(npix))

! import input map
call input_map(fin, Min, npix, nmaps)

! import mask if specified
if (nArguments() < 3) then; Mask = 1.0; nuse = npix; else
	call getArgument(3, fmask)
	call input_map(fmask, Mask, npix, nmaps)
	
	! masked pixels are not ranked
	nuse = 0; do i = 1,npix
		if (Mask(i,ch) == 0.0) then
			Min(i,ch) = HUGE(Min)
		else
			nuse = nuse + 1
		end if
	end do
end if

! compute rank ordering and L-weights
call gegenbauer(P, nmoms-1, 1.0)
call indexx(npix, Min(:,ch), indx)
call rankx(npix, indx, rank)

do i = 1,npix
        Mout(i,:) = Mask(i,ch) * matmul(P, X(rank(i), nuse, nmoms-1))
end do

! output L-weight masks (to a single FITS container)
call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator='LMASK', version='$Revision$')
call output_map(Mout, header, '!'//fout)

! output L-weighted maps (to separate FITS files) if requested
if (nArguments() > 3) then
	call getArgument(4, fout)
	
	do i = 1,npix
		Mout(i,:) = Min(i,ch) * Mout(i,:)
	end do
	
	do i = 1,nmoms
		write (fmask,'(A,A1,I1,A5)') trim(fout), '-', i, '.fits'
		
		call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator='LMASK', version='$Revision$')
		call output_map(Mout(:,(/i/)), header, '!'//fmask)
	end do
end if


contains

! scaled Gegenbauer polynomials P(k,alpha,x) = C(k,alpha/2,x)/C(k,alpha/2,1)
! generalize Legendre P (alpha=1), Chebyshev T (alpha=0) and Chebyshev U (alpha=2)
! satisfy P(0) = 1.0; P(1) = x; P(k+1) = ((2*k+alpha)*x*P(k) - k*P(k-1))/(k+alpha)
! returns *shifted* polynomial coefficients as P(k,alpha,2*x-1) = sum(P(k,n)*x^n)
subroutine gegenbauer(P, l, alpha)
        integer k, l; real P(0:l,0:l), alpha
        
        P = 0.0; P(0,0) = 1.0; P(1,1) = 2.0; P(1,0) = -1.0; do k = 1,l-1
                P(k+1,:) = ((2*k+alpha)*(2.0*cshift(P(k,:),-1)-P(k,:)) - k*P(k-1,:))/(k+alpha)
        end do
end subroutine gegenbauer

! Calculate L-weights from rank order vector X using matmul(P,X)
function X(r,n,l)
	integer r, n, l, k; real X(0:l)
	
	X(0) = 1.0; do k = 1,l; X(k) = X(k-1) * (r-k)/(n-k); end do
end function X

end
