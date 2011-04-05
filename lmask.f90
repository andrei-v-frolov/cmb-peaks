! $Id$
! Synthesize a rank-ordered weigth masks for calculating L-moments
! invoke: lmask <map.fits[:channel]> <lmask-base[:moments]> [mask.fits]

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
real(DP), allocatable :: P(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! input map and optional signal channel selection
call getArgument(1, fin ); i = index(fin,  ":", .true.)
if (i > 0) then; read (fin(i+1:),*) ch; fin(i:) = ""; end if

! output map and optional number of moments selection
call getArgument(2, fout); i = index(fout, ":", .true.)
if (i > 0) then; read (fout(i+1:),*) nmoms; fout(i:) = ""; end if

! allocate dynamic arrays to store maps and ranks
ntot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord); npix = nside2npix(nside)
allocate(Min(npix,nmaps), Mask(npix,nmaps), Mout(npix,nmoms), P(0:nmoms), indx(npix), rank(npix))

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
call indexx(npix, Min(:,ch), indx)
call rankx(npix, indx, rank)

do i = 1,npix
        call legendre(P, rank(i), nuse, nmoms-1)
        Mout(i,:) = Mask(i,ch) * P(1:nmoms)
end do

! output L-weight masks (to a single FITS container)
!call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator='LMASK', version='$Revision$')
!call output_map(Mout, header, '!'//fout)

! output L-weight masks (to separate FITS files)
do i = 1,nmoms
	write (fmask,'(A,A1,I1,A5)') trim(fout), '-', i, '.fits'
	
	call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator='LMASK', version='$Revision$')
	call output_map(Mout(:,(/i/)), header, '!'//fmask)
end do

contains

! Calculate Legendre weights from rank ordering
subroutine legendre(P, r, n, l)
        integer k, r, n, l; real(DP) x, P(-1:l)
        
        P(-1:0) = 1.0; do k = 1,l; x = real(r-k)/real(n-k)
                P(k) = ((2*k-1)*(2.0*x-1.0)*P(k-1) - (k-1)*P(k-2))/k
        end do
end subroutine legendre

end
