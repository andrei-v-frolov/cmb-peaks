! $Id$
! Synthesize a rank-ordered weigth masks for calculating L-moments
! invoke: lmask <map.fits[:channel]> <lmask.fits[:moments]> [mask.fits[:channel]] [lmap-base]

program lmask

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: datach = 1, maskch = 1, nmoms = 4    ! defaults
integer :: i, npix, nuse, nside = 0, ord = 0    ! map format

character(len=80) :: header(64), fin, fout, fmask
real(SP), allocatable :: Min(:), Mask(:), Mout(:,:)
integer, allocatable :: indx(:), rank(:)
real(DP), allocatable :: P(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! input map (and optional signal channel selection)
! output map (and optional number of moments selection)
call getArgument(1, fin ); call parse(fin, datach)
call getArgument(2, fout); call parse(fout, nmoms)

! read input map
call read_channel(fin, Min, nside, datach, ord); npix = nside2npix(nside)

! read mask if specified
if (nArguments() < 3) then
        allocate(Mask, mold=Min)
        Mask = 1.0; nuse = npix
else
        call getArgument(3, fmask); call parse(fmask, maskch)
	call read_channel(fmask, Mask, nside, maskch, ord)
	
	! masked pixels are not ranked
	nuse = 0; do i = 1,npix
		if (Mask(i) == 0.0) then
			Min(i) = HUGE(Min)
		else
			nuse = nuse + 1
		end if
	end do
end if

! allocate dynamic arrays to store output maps and ranks
allocate(Mout(npix,nmoms), P(nmoms,nmoms), indx(npix), rank(npix))

! compute rank ordering and L-weights
call gegenbauer(P, nmoms-1, 1.0)
call indexx(npix, Min, indx)
call rankx(npix, indx, rank)

do i = 1,npix
        Mout(i,:) = Mask(i) * matmul(P, X(rank(i), nuse, nmoms-1))
end do

! output L-weight masks (to a single FITS container)
call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator='LMASK', version='$Revision$')
call output_map(Mout, header, '!'//fout)

! output L-weighted maps (to separate FITS files) if requested
if (nArguments() > 3) then
	call getArgument(4, fout)
	
	do i = 1,npix
		Mout(i,:) = Min(i) * Mout(i,:)
	end do
	
	do i = 1,nmoms
		write (fmask,'(A,A1,I1,A5)') trim(fout), '-', i, '.fits'
		
		call write_minimal_header(header, 'MAP', nside=nside, order=ord, creator='LMASK', version='$Revision$')
		call output_map(Mout(:,(/i/)), header, '!'//fmask)
	end do
end if


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse optional channel specification
subroutine parse(file, channel)
        character(*) file
        integer channel, i
        
        i = index(file,  ":", .true.)
        
        if (i > 0) then
                read (file(i+1:),*) channel
                file(i:) = ""
        end if
end subroutine parse

! read a single map channel, allocating storage if necessary
subroutine read_channel(fin, M, nside, channel, ord)
        character(*) fin
        real(SP), allocatable :: M(:), TMP(:,:)
        integer channel, nside, npix, nmaps, ord
        
        ! read full map into temporary storage
        nmaps = 0; call read_map(fin, TMP, nside, nmaps, ord)
        if (channel > nmaps) call abort(trim(fin) // ": too few channels in an input map")
        
        ! allocate storage if needed
        npix = nside2npix(nside); if (.not. allocated(M)) allocate(M(npix))
        if (size(M) /= npix) call abort(trim(fin) // ": unexpected storage array shape")
        
        ! copy over the data we want, free the full map
        M = TMP(:,channel); deallocate(TMP)
end subroutine read_channel

! read map from FITS file, allocating storage if necessary
subroutine read_map(fin, M, nside, nmaps, ord)
        character(*) fin
        real(SP), allocatable :: M(:,:)
        integer nside, npix, nmaps, ord
        
        ! read header info
        character(len=80) :: header(64)
        integer hside, htot, hmaps, hord
        integer, parameter :: RING = 1, NEST = 2
        
        htot = getsize_fits(fin, nside=hside, nmaps=hmaps, ordering=hord)
        if (htot == -1) call abort(trim(fin) // ": file not found")
        
        ! check if map format agrees with requested one
        if (nside == 0) nside = hside; if (hside /= nside) call abort(trim(fin) // ": map resolution does not conform")
        if (nmaps == 0) nmaps = hmaps; if (hmaps  < nmaps) call abort(trim(fin) // ": too few channels in an input map")
        if (  ord == 0)   ord = hord;  if ( hord /= ord) call warning(trim(fin) // ": map ordering is being converted")
        
        ! allocate storage if needed
        npix = nside2npix(nside); if (.not. allocated(M)) allocate(M(npix,nmaps))
        if (size(M,1) /= npix .or. size(M,2) < nmaps) call abort(trim(fin) // ": unexpected storage array shape")
        
        ! read map data, converting order if needed
        call input_map(fin, M, npix, nmaps)
        if (hord == RING .and. ord == NEST) call convert_ring2nest(nside, M)
        if (hord == NEST .and. ord == RING) call convert_nest2ring(nside, M)
end subroutine read_map

! output warning message to stderr
subroutine warning(msg)
        character(*) msg; logical, parameter :: verbose = .true.
        
        if (verbose) write (0,'(a)') "warning: " // msg
end subroutine warning

end
