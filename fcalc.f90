! $Id$
! HEALPix map calculator, produces output map from zero or more inputs via
!  prefix operator: fcalc 'x' M.fits [=:] output.fits
! postfix operator: fcalc M.fits 'x' [=:] output.fits
!  binary operator: fcalc M1.fits 'x' M2.fits [=:] output.fits
! ternary operator: fcalc M1.fits 'x' M2.fits 'y' M3.fits [=:] output.fits
!  correlated maps: fcalc 'x' M1.fits x M2.fits [=:] out1.fits x out2.fits
! see source code for complete list of operators currently implemented

program fcalc

! HEALPix includes
use mapio
use pdetools
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

character(len=8000) :: op, fin1, fin2, fin3, fout, fout2
integer :: nmaps = 0, nside = 0, ord = 0, pol = -1
integer i, j, n, lmin, lmax, bands, seed(2)

real(IO), dimension(:,:), allocatable :: M1, M2, M3, Mout, Mout2
logical, dimension(:,:), allocatable :: valid
real, allocatable :: bandpass(:,:)
real multipoles(0:3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
if (.not. (prefix() .or. postfix() .or. binary() .or. ternary() .or. split())) call abort("cannot parse command line expression supplied")
if (.not. allocated(M1)) call abort("no map data supplied, I am done")

! map parameters
n = nside2npix(nside)-1; lmax = 3*nside-1; lmin = 1

! output storage
allocate(Mout, mold=M1); Mout = 0.0

! valid data mask
allocate(valid(0:n,nmaps), source=.true.)
if (allocated(M1) .and. .not. allocated(M2) .and. .not. allocated(M3)) valid = .not. (isnan(M1))
if (allocated(M1) .and. allocated(M2) .and. .not. allocated(M3)) valid = .not. (isnan(M1) .or. isnan(M2))
if (allocated(M1) .and. allocated(M2) .and. allocated(M3)) valid = .not. (isnan(M1) .or. isnan(M2) .or. isnan(M3))

! initialize random number generator (use urandom on clusters!)
if (index(op,'randomize') > 0 .or. index(op,'shuffle') > 0) then
	open (333, file="/dev/random", action='read', form='binary')
	read (333) seed; call random_seed(PUT=seed); close (333)
	!seed = [123456789,987654321]; call random_seed(PUT=seed)
end if

! apply operator
select case (op)
	! arithmetic operators
	case ('+');  Mout = M1 + M2
	case ('-');  Mout = M1 - M2
	case ('*');  Mout = M1 * M2
	case ('/');  Mout = M1 / M2
	case ('//'); Mout = floor(M1/M2)
	case ('**'); Mout = M1 ** M2
	case ('sqrt'); Mout = sqrt(M1)
	case ('accumulate'); Mout = M1 + M2*M2
	case ('accumulate-'); Mout = M1 + (M2-M3)**2
	
	! comparison operators
	case ('<');  where (M1 < M2) Mout = 1.0
	case ('>');  where (M1 > M2) Mout = 1.0
	case ('<='); where (M1 <= M2) Mout = 1.0
	case ('>='); where (M1 >= M2) Mout = 1.0
	case ('=','=='); where (M1 == M2) Mout = 1.0
	case ('!=','/=','<>'); where (M1 /= M2) Mout = 1.0
	
	! rank-order map, outputing CDF value for valid pixels (on per channel basis)
	case ('rank'); do i = 1,nmaps; call percentile(nside, M1(:,i), valid(:,i), Mout(:,i)); end do
	
	! projection operators
	case ('project on'); Mout = sum(M1*M2,valid)/sum(M2*M2,valid) * M2
	case ('orthogonal'); Mout = M1 - sum(M1*M2,valid)/sum(M2*M2,valid) * M2
	case ('remove monopole'); if (pol > 0) call abort(trim(op) // " operator works on scalar maps only")
		do i = 1,nmaps; call remove_dipole(nside, M1(:,i), ord, 1, multipoles, [-1.0,1.0], mask=M2(:,i)); end do; Mout = M1
	case ('remove dipole'); if (pol > 0) call abort(trim(op) // " operator works on scalar maps only")
		do i = 1,nmaps; call remove_dipole(nside, M1(:,i), ord, 2, multipoles, [-1.0,1.0], mask=M2(:,i)); end do; Mout = M1
	
	! projection operators (masked version)
	case ('project onwith'); Mout = sum(M1*M2*M3,valid)/sum(M2*M2*M3,valid) * M2
	case ('orthogonalwith'); Mout = M1 - sum(M1*M2*M3,valid)/sum(M2*M2*M3,valid) * M2
	
	! masking operators
	case ('valid'); where (valid) Mout = 1.0
	case ('invalid'); where (.not. valid) Mout = 1.0
	case ('mask'); Mout = M1*M2; where (M2 == 0.0) Mout = 1.0/0.0
	case ('unmask'); Mout = M1/M2; where (M2 == 0.0) Mout = 1.0/0.0
	case ('within:'); where (M1 >= M2 .and. M1 <= M3) Mout = 1.0
	case ('apodize:'); Mout = apodize((M1-M2)/(M3-M2))
	
	! mask manipulation
	case ('grow'); do i = 1,nmaps; call grow_mask(nside, ord, M1(:,i), Mout(:,i)); end do
	case ('shrink'); do i = 1,nmaps; call shrink_mask(nside, ord, M1(:,i), Mout(:,i)); end do
	case ('sources'); do i = 1,nmaps; call sources_mask(nside, ord, 3.0, 60.0, M1(:,i), Mout(:,i)); end do
	
	! inpainting and filling
	case ('inpaint'); do i = 1,nmaps; call inpaint(nside, ord, M1(:,i), M2(:,i), Mout(:,i)); end do
	case ('inpaintmass'); do i = 1,nmaps; call inpaint(nside, ord, M1(:,i), M2(:,i), Mout(:,i), m2=M3(:,i)); end do
	case ('inpaintwith'); do i = 1,nmaps; call inpaint(nside, ord, M1(:,i), M2(:,i), Mout(:,i), fill=M3(:,i)); end do
	case ('inpaintapodize'); do i = 1,nmaps; call inpaint(nside, ord, M1(:,i), M2(:,i), Mout(:,i), apo=M3(:,i)); end do
	
	! conversion operators
	case ('nest');
		select case (ord)
			case (NEST,0); Mout = M1; ord = NEST
			case (RING);   Mout = M1; call convert_ring2nest(nside, Mout); ord = NEST
			case default; call abort(trim(op) // " conversion encountered unkown ordering")
		end select
	case ('ring');
		select case (ord)
			case (RING,0); Mout = M1; ord = RING
			case (NEST);   Mout = M1; call convert_nest2ring(nside, Mout); ord = RING
			case default; call abort(trim(op) // " conversion encountered unkown ordering")
		end select
	
	! reduction operators
	case ('any');     nmaps = 1; pol = -1; where (any(M1 /= 0.0,2)) Mout(:,1) = 1.0
	case ('all');     nmaps = 1; pol = -1; where (all(M1 /= 0.0,2)) Mout(:,1) = 1.0
	case ('sum');     nmaps = 1; pol = -1; Mout(:,1) = sum(M1,2)
	case ('norm');    nmaps = 1; pol = -1; Mout(:,1) = sqrt(sum(M1**2,2))
	case ('product'); nmaps = 1; pol = -1; Mout(:,1) = product(M1,2)
	case ('select');  nmaps = 1; pol = -1; forall (i=0:n) Mout(i,1) = M1(i,M2(i,1))
	
	! composition operators
	case ('zip');     if (nmaps /= 1) call abort(trim(op) // " composition works on single channel maps")
		nmaps=2; deallocate(Mout); allocate(Mout(0:n,nmaps)); Mout(:,1) = M1(:,1); Mout(:,2) = M2(:,1)
	case ('zipwith'); if (nmaps /= 1) call abort(trim(op) // " composition works on single channel maps")
		nmaps=3; deallocate(Mout); allocate(Mout(0:n,nmaps)); Mout(:,1) = M1(:,1); Mout(:,2) = M2(:,1); Mout(:,3) = M3(:,1)
	
	! injection operators
	case ('replacewith'); Mout = M1; forall (i=0:n) Mout(i,M2(i,1)) = M3(i,M2(i,1))
	case ('updatewith');  Mout = M1; forall (i=0:n) Mout(i,M2(i,1)) = M1(i,M2(i,1)) + M3(i,M2(i,1))
	case ('update+'); Mout = M1; forall (i=0:n) Mout(i,M2(i,1)) = M1(i,M2(i,1)) + M3(i,M2(i,1))
	case ('update-'); Mout = M1; forall (i=0:n) Mout(i,M2(i,1)) = M1(i,M2(i,1)) - M3(i,M2(i,1))
	case ('update*'); Mout = M1; forall (i=0:n) Mout(i,M2(i,1)) = M1(i,M2(i,1)) * M3(i,M2(i,1))
	case ('update/'); Mout = M1; forall (i=0:n) Mout(i,M2(i,1)) = M1(i,M2(i,1)) / M3(i,M2(i,1))
	
	! randomize operators
	case ('randomize');
		allocate(M2, mold=M1); allocate(M3, mold=M1)
		call random_number(M2); call random_number(M3)
		Mout = M1 * sqrt(-2.0*log(M2)) * cos(2.0*pi*M3)
	case ('shuffle');
		allocate(M2, mold=M1); call random_number(M2)
		do i = 0,n; j = n - floor((n-i+1)*M2(i,1))
			Mout(i,:) = M1(j,:); M1(j,:)=M1(i,:)
		end do
	case ('randomize-alm');
		select case (nmaps)
			case (1); call randomize(nside, ord, lmin, lmax, M1(:,1), Mout(:,1))
			case (3); call randomize(nside, ord, lmin, lmax, M1(:,1), Mout(:,1))
			          call randomize_qu(nside, ord, lmin, lmax, M1(:,2:3), Mout(:,2:3), .true., .true.)
			case (2); call randomize_qu(nside, ord, lmin, lmax, M1(:,1:2), Mout(:,1:2), .true., .true.)
			case default; call abort(trim(op) // " conversion requires I, QU or IQU map format")
		end select
	case ('randomize-blm');
		select case (nmaps)
			case (2); call randomize_qu(nside, ord, lmin, lmax, M1(:,1:2), Mout(:,1:2), randomizeB=.true.)
			case (3); call randomize_qu(nside, ord, lmin, lmax, M1(:,2:3), Mout(:,2:3), randomizeB=.true.); Mout(:,1) = M1(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	case ('randomize-elm');
		select case (nmaps)
			case (2); call randomize_qu(nside, ord, lmin, lmax, M1(:,1:2), Mout(:,1:2), randomizeE=.true.)
			case (3); call randomize_qu(nside, ord, lmin, lmax, M1(:,2:3), Mout(:,2:3), randomizeE=.true.); Mout(:,1) = M1(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	
	! cross-correlated randomize
	case ('xrandomize-blm');
		select case (nmaps)
			case (2); call xrandomize_qu(nside, ord, lmin, lmax, M1(:,1:2), M2(:,1:2), Mout(:,1:2), Mout2(:,1:2), randomizeB=.true.)
			case (3); call xrandomize_qu(nside, ord, lmin, lmax, M1(:,2:3), M2(:,2:3), Mout(:,2:3), Mout2(:,2:3), randomizeB=.true.)
				  Mout(:,1) = M1(:,1); Mout2(:,1) = M2(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	case ('xrandomize-elm');
		select case (nmaps)
			case (2); call xrandomize_qu(nside, ord, lmin, lmax, M1(:,1:2), M2(:,1:2), Mout(:,1:2), Mout2(:,1:2), randomizeE=.true.)
			case (3); call xrandomize_qu(nside, ord, lmin, lmax, M1(:,2:3), M2(:,2:3), Mout(:,2:3), Mout2(:,2:3), randomizeE=.true.)
				  Mout(:,1) = M1(:,1); Mout2(:,1) = M2(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	
	! one-point operators
	case ('frac');
		select case (nmaps)
			case (3); forall (i=0:n) Mout(i,:) = fraction(M1(i,:))
			case default; call abort(trim(op) // " conversion requires IQU map format")
		end select
	case ('log');
		select case (nmaps)
			case (1); Mout = log(M1)
			case (3); forall (i=0:n) Mout(i,:) = log_iqu(M1(i,:))
			case default; call abort(trim(op) // " conversion requires I or IQU map format")
		end select
	case ('exp');
		select case (nmaps)
			case (1); Mout = exp(M1)
			case (3); forall (i=0:n) Mout(i,:) = exp_iqu(M1(i,:))
			case default; call abort(trim(op) // " conversion requires I or IQU map format")
		end select
	
	! vector field operators
	case ('cartesian->healpix','xyz->XYZ');
		select case (nmaps)
			case (3); do i = 0,n; Mout(i,:) = cart2hlpx(nside, ord, i, M1(i,:), +1); end do
			case default; call abort(trim(op) // " reconstruction requires B[xyz] map as input")
		end select
	case ('healpix->cartesian','XYZ->xyz');
		select case (nmaps)
			case (3); do i = 0,n; Mout(i,:) = cart2hlpx(nside, ord, i, M1(i,:), -1); end do
			case default; call abort(trim(op) // " reconstruction requires B[XYZ] map as input")
		end select
	case ('XY->EB');
		select case (nmaps)
			case (2); call rotate_qu2eb(nside, ord, lmax, M1(:,1:2), Mout(:,1:2), spin=1)
			case (3); call rotate_qu2eb(nside, ord, lmax, M1(:,1:2), Mout(:,1:2), spin=1); Mout(:,3) = M1(:,3)
			case default; call abort(trim(op) // " conversion requires XY or XYZ map format")
		end select
	case ('EB->XY');
		select case (nmaps)
			case (2); call rotate_eb2qu(nside, ord, lmax, M1(:,1:2), Mout(:,1:2), spin=1)
			case (3); call rotate_eb2qu(nside, ord, lmax, M1(:,1:2), Mout(:,1:2), spin=1); Mout(:,3) = M1(:,3)
			case default; call abort(trim(op) // " conversion requires EB or EBZ map format")
		end select
	
	! polarization operators
	case ('QU->EB');
		select case (nmaps)
			case (2); call rotate_qu2eb(nside, ord, lmax, M1(:,1:2), Mout(:,1:2))
			case (3); call rotate_qu2eb(nside, ord, lmax, M1(:,2:3), Mout(:,2:3)); Mout(:,1) = M1(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	case ('EB->QU');
		select case (nmaps)
			case (2); call rotate_eb2qu(nside, ord, lmax, M1(:,1:2), Mout(:,1:2))
			case (3); call rotate_eb2qu(nside, ord, lmax, M1(:,2:3), Mout(:,2:3)); Mout(:,1) = M1(:,1)
			case default; call abort(trim(op) // " conversion requires EB or IEB map format")
		end select
	case ('QU->pure EB');
		select case (nmaps)
			case (2); call rotate_qu2eb_pure(nside, ord, lmax, M1(:,1:2), M2(:,1), Mout(:,1:2))
			case (3); call rotate_qu2eb_pure(nside, ord, lmax, M1(:,2:3), M2(:,2), Mout(:,2:3)); Mout(:,1) = M1(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	case ('inpaint QU');
		select case (nmaps)
			case (2); call inpaint_qu(nside, ord, M1(:,1:2), M2(:,1), Mout(:,1:2))
			case (3); call inpaint_qu(nside, ord, M1(:,2:3), M2(:,2), Mout(:,2:3))
			          call inpaint(nside, ord, M1(:,1), M2(:,1), Mout(:,1))
			case default; call abort(trim(op) // " tensor inpainting requires QU or IQU map format")
		end select
	case ('purify');
		select case (nmaps)
			case (2); call inpaint_purified_qu(nside, ord, lmax, M1(:,1:2), M2(:,1), Mout(:,1:2))
			case (3); call inpaint_purified_qu(nside, ord, lmax, M1(:,2:3), M2(:,2), Mout(:,2:3))
			          call inpaint(nside, ord, M1(:,1), M2(:,1), Mout(:,1))
			case default; call abort(trim(op) // " purified inpainting requires QU or IQU map format")
		end select
	
	! spectral operators
	case ('bbodyfrequency');
		if (nmaps /=1) call warning("modified blackbody correction works on single channel maps")
		nmaps = 1; forall (i=0:n) Mout(i,1) = bbody_log_tcmb(real(M3(i,1)), real(M1(i,1)), real(M2(i,1)))
	
	case ('bbodybandpass');
		call read_bandpass(bandpass, bands)
		if (nmaps /=1) call warning("modified blackbody correction works on single channel maps")
		if (verbose) write (*,*) "Computing modified blackbody color correction, sample output:"
		nmaps = 1; do i = 0,n
			Mout(i,1) = bbody_log_cc(real(M3(i,1)), real(M1(i,1)), real(M2(i,1)), bandpass, bands)
			if (verbose .and. mod(i,1048576) == 0) write (*,*) i, Mout(i,1)
		end do
	
	! reconstruction operators
	case ('magnetic->pqu');
		select case (nmaps)
			case (3); do i = 0,n; Mout(i,:) = magnetic2pqu(nside, ord, i, M1(i,:)); end do
			case default; call abort(trim(op) // " reconstruction requires B[xyz] map as input")
		end select
	case ('pqu->magnetic');
		select case (nmaps)
			case (3); call magnetic_fit(nside, ord, 0, 0, M1(:,2:3), Mout)
			case default; call abort(trim(op) // " reconstruction requires pqu map as input")
		end select
	case ('pqu->magneticlmax');
		select case (nmaps)
			case (3); call magnetic_fit(nside, ord, 0, nint(M3(0,1)), M1(:,2:3), Mout, M2)
			case default; call abort(trim(op) // " reconstruction requires pqu map as input")
		end select
	
	! unknown operator
	case default; call abort(trim(op) // ": operation not supported")
end select

! write output map(s)
call write_map(fout, Mout(:,1:nmaps), nside, ord, pol, creator='FCALC')
if (allocated(Mout2)) call write_map(fout2, Mout2(:,1:nmaps), nside, ord, pol, creator='FCALC')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse prefix operator command line
function prefix()
	character(len=80) :: x; logical prefix; prefix = .false.
	
	! argument count guard
	if (nArguments() < 3 .or. nArguments() > 4) return
	
	! operator placement
	call getArgument(1, x)
	
	! prefix operation guard
	select case (x)
		case ('frac','log','exp','rank','sqrt','valid','invalid','any','all','sum','norm','product')
		case ('randomize','shuffle','randomize-alm','randomize-blm','randomize-elm')
		case ('xyz->XYZ','XYZ->xyz','XY->EB','EB->XY','QU->EB','EB->QU')
		case default; return
	end select
	
	! operator name
	prefix = .true.; op = trim(x)
	
	! read input maps
	call getArgument(2, fin1); call read_map(fin1, M1, nside, nmaps, ord, pol)
	
	! output map name
	call getArgument(3, fout); if (fout .eq. '=:') call getArgument(4, fout)
end function

! parse postfix operator command line
function postfix()
	character(len=80) :: x; logical postfix; postfix = .false.
	
	! argument count guard
	if (nArguments() < 3 .or. nArguments() > 4) return
	
	! operator placement
	call getArgument(2, x)
	
	! postfix operation guard
	select case (x)
		case ('nest','ring','grow','shrink','sources','XY->EB','EB->XY','QU->EB','EB->QU')
		case ('cartesian->healpix','healpix->cartesian','magnetic->pqu','pqu->magnetic')
		case default; return
	end select
	
	! operator name
	postfix = .true.; op = trim(x)
	
	! read input maps
	call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord, pol)
	
	! output map name
	call getArgument(3, fout); if (fout .eq. '=:') call getArgument(4, fout)
end function

! parse binary operator command line
function binary()
	character(len=80) :: x; logical binary; binary = .false.
	
	! argument count guard
	if (nArguments() < 4 .or. nArguments() > 5) return
	
	! operator placement
	call getArgument(2, x)
	
	! binary operation guard
	select case (x)
		case ('+','-','*','/','//','**')
		case ('<','>','<=','>=','=','==','!=','/=','<>')
		case ('project on','orthogonal','remove monopole','remove dipole','accumulate','select','zip')
		case ('valid','invalid','mask','unmask','inpaint','inpaint QU','QU->pure EB','purify')
		case default; return
	end select
	
	! operator name
	binary = .true.; op = trim(x)
	
	! read input maps
	call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord, pol)
	call getArgument(3, fin2); call read_map(fin2, M2, nside, nmaps, ord, pol)
	
	! output map name
	call getArgument(4, fout); if (fout .eq. '=:') call getArgument(5, fout)
end function

! parse ternary operator command line
function ternary()
	character(len=80) :: x, y; logical ternary; ternary = .false.
	
	! argument count guard
	if (nArguments() < 6 .or. nArguments() > 7) return
	
	! operator placements
	call getArgument(2, x)
	call getArgument(4, y)
	
	! is it really ternary?
	select case (x)
		case ('zip','replace'); if (y .eq. 'with') ternary = .true.
		case ('update');  select case (y); case ('with','+','-','*','/'); ternary = .true.; end select
		case ('inpaint'); select case (y); case ('mass','with','apodize'); ternary = .true.; end select
		case ('project on','orthogonal'); if (y .eq. 'with') ternary = .true.
		case ('accumulate'); if (y .eq. '-') ternary = .true.
		case ('within','apodize'); if (y .eq. ':') ternary = .true.
		case ('pqu->magnetic'); if (y .eq. 'lmax') ternary = .true.
		case ('bbody'); if (y .eq. 'frequency' .or. y .eq. 'bandpass') ternary = .true.
	end select
	
	! ternary operation guard
	if (.not. ternary) return; op = trim(x) // trim(y)
	
	! read input maps
	call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord, pol)
	call getArgument(3, fin2); call read_map(fin2, M2, nside, nmaps, ord, pol)
	call getArgument(5, fin3); call read_map(fin3, M3, nside, nmaps, ord, pol)
	
	! output map name
	call getArgument(6, fout); if (fout .eq. '=:') call getArgument(7, fout)
end function

! parse 2x2 split operator command line
function split()
	character(len=80) :: x, y, z; logical split; split = .false.
	
	! argument count guard
	if (nArguments() < 7 .or. nArguments() > 8) return
	
	! operator placements
	call getArgument(1, x)
	call getArgument(3, y)
	call getArgument(nArguments()-1, z)
	if (y /= 'x' .or. z /= 'x') return
	
	! is it really split?
	select case (x)
		case ('randomize-alm','randomize-blm','randomize-elm')
		case default; return
	end select
	
	! operator name
	split = .true.; op = 'x' // trim(x)
	
	! read input maps
	call getArgument(2, fin1); call read_map(fin1, M1, nside, nmaps, ord, pol)
	call getArgument(4, fin2); call read_map(fin2, M2, nside, nmaps, ord, pol)
	
	! output map names
	call getArgument(5, fout)
	if (fout .eq. '=:') then
		call getArgument(6, fout)
		call getArgument(8, fout2)
	else
		call getArgument(7, fout2)
	end if
	
	! allocate second output map
	allocate(Mout2, mold=M1); Mout2 = 0.0
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! apodization function smoothly interpolates between 0 and 1
elemental function apodize(x)
	real(IO) x, apodize; intent(in) x
	
	apodize = 0.0; if (x <= 0.0) return
	apodize = 1.0; if (x >= 1.0) return
	
	!apodize = (1.0+tanh(tan(pi*(x-0.5))))/2.0
	apodize = sin(pi/2.0*x)**2
end function

! polarization fraction
pure function fraction(iqu)
	real(IO) fraction(3), iqu(3); intent(in) iqu
	
	associate (I => real(iqu(1),DP), Q => real(iqu(2),DP), U => real(iqu(3),DP))
		fraction = iqu/(I + sqrt(Q*Q+U*U))
	end associate
end function

! logarithm of polarization tensor
pure function log_iqu(iqu)
	real(IO) log_iqu(3), iqu(3); intent(in) iqu
	real(DP) P, PxP
	
	associate (I => real(iqu(1),DP), Q => real(iqu(2),DP), U => real(iqu(3),DP))
		PxP = Q*Q + U*U; P = sqrt(PxP)
		log_iqu(1) = log(I*I - PxP)/2.0
		log_iqu(2:3) = [Q,U]/P * log((I+P)/(I-P))/2.0
	end associate
end function

! exponent of polarization tensor
pure function exp_iqu(iqu)
	real(IO) exp_iqu(3), iqu(3); intent(in) iqu
	real(DP) P, PxP
	
	associate (I => real(iqu(1),DP), Q => real(iqu(2),DP), U => real(iqu(3),DP))
		PxP = Q*Q + U*U; P = sqrt(PxP)
		exp_iqu(1) = exp(I) * cosh(P)
		exp_iqu(2:3) = [Q,U]/P * exp(I) * sinh(P)
	end associate
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! rank-order map, outputing CDF value for valid pixels
subroutine percentile(nside, map, valid, cdf)
	integer nside, npix, used
	real(IO), dimension(0:12*nside**2-1) :: map, cdf
	logical, dimension(0:12*nside**2-1) :: valid
	real(DP), allocatable :: M(:)
	integer, allocatable :: idx(:), rank(:)
	
	npix = nside2npix(nside)
	allocate(M(npix), idx(npix), rank(npix))
	
	where (.not. valid) map = HUGE(map)
	M = map; used = count(valid)
	
	call indexx(npix, M, idx)
	call rankx(npix, idx, rank)
	
	cdf = (rank-1)/(used-1.0)
	where (.not. valid) cdf = 1.0/0.0
	
	deallocate(M, idx, rank)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! hyperbolic sinc function
elemental function sinch(x)
	real sinch, x; intent(in) x
	
	sinch = 1.0; if (x /= 0.0) sinch = sinh(x)/x
end function

! logarithm of hyperbolic sinc function, accurate for all argument values
elemental function lsinch(x)
	real lsinch, x, p, q; intent(in) x
	
	real, parameter :: c(7) = [ &
		0.4132452269851979093654156131033060008760029088251892994339693444Q-1, &
		0.8478301007195085738735669818781708187733967110608704491525621626Q-4, &
		3.319067386208102157904215623529367549226693561719962519708929754Q-07, &
		1.535892374193619383767933303677575365217042971636305586372502646Q-09, &
		7.659483473137847811201375808541567436659033038955381198516298508Q-12, &
		3.988505147571547603481766897703709232256815977906297345067145370Q-14, &
		2.137519596232463056459672303515036468045220079963504661934751418Q-16  &
	]
	
	p = 2.0*abs(x); q = exp(-p)
	
	if (p < 1.0) then
		lsinch = sum(c * sin([1:7]*asin(p))**2)
	else
		lsinch = p/2.0 - log(p) - 2.0*atanh(q/(2.0-q))
	end if
end function

! modified blackbody spectrum, in T[RJ] units
pure function bbody_log_trj(nu, t, beta, nu0)
	real bbody_log_trj, nu, mu, t, beta, gamma, nu0
	intent (in) nu, t, beta, nu0; optional nu0
	
	! h/2kT, in units of [1/GHz]
	gamma = 0.023996223311/t
	
	! reference frequency, 353GHz by default
	mu = 353.0; if (present(nu0)) mu = nu0
	
	bbody_log_trj = beta*log(nu/mu) + gamma*(mu-nu) + (lsinch(gamma*mu) - lsinch(gamma*nu))
end function

! modified blackbody spectrum, in T[CMB] units
pure function bbody_log_tcmb(nu, t, beta, nu0)
	real bbody_log_tcmb, nu, mu, t, beta, nu0
	intent(in) nu, t, beta, nu0; optional nu0
	
	! COBE/FIRAS CMB temperature is 2.7255K
	real, parameter :: gamma = 0.023996223311/2.7255
	
	! reference frequency, 353GHz by default
	mu = 353.0; if (present(nu0)) mu = nu0
	
	bbody_log_tcmb = 2.0*(lsinch(gamma*nu) - lsinch(gamma*mu)) + bbody_log_trj(nu, t, beta, mu)
end function

! read bandpass correction data
subroutine read_bandpass(b, n)
	real, allocatable :: b(:,:); integer i, n, status
	real, parameter :: c = 29.9792458 ! speed of light, cm*GHz
	
	! allocate bandpass storage if not already done
	if (allocated(b)) then; n = size(b,2); else; n = 2001; allocate(b(2,n)); end if
	
	! parse bandpass data (wavenumber [1/cm], transmission [n/a])
	i = 1; do while (i <= n)
		read (*,*,iostat=status) b(:,i)
		if (status < 0) exit
		if (status > 0) cycle
		i = i + 1
	end do
	
	! sanity checks on bandpass data read
	if (i > n) call warning("bandpass data truncated to fit allocated storage!")
	if (i < 2) call abort("no bandpass data supplied!"); n = i - 1
	
	! convert frequency units to GHz
	b(1,:) = b(1,:) * c
end subroutine

! simple trapezoid rule integrator for tabulated functions
function integral(v, x, n)
	real integral, v(n), x(n); integer i, n; real(DP) s
	
	s = 0.0; do i = 1,n-1
		s = s + (v(i+1) + v(i)) * (x(i+1) - x(i))
	end do
	
	integral = s/2.0
end function

! bandpass integral of dB_nu(T)/dT referred to central frequency mu
function bandpass_dBdT(mu, t, beta, bandpass, n)
	real bandpass_dBdT, mu, t, beta, gamma, bandpass(:,:); integer n
	
	! h/2kT, in units of [1/GHz]
	gamma = 0.023996223311/t
	
	associate(nu => bandpass(1,1:n), tau => bandpass(2,1:n))
	bandpass_dBdT = integral(tau * nu**(beta+2.0)/sinch(gamma*nu)**2, nu, n) * sinch(gamma*mu)**2/mu**(beta+2.0)
	end associate
end function

! bandpass integral of B_nu(T) referred to central frequency mu
function bandpass_dBdA(mu, t, beta, bandpass, n)
	real bandpass_dBdA, mu, t, beta, gamma, bandpass(:,:); integer n
	
	! h/2kT, in units of [1/GHz]
	gamma = 0.023996223311/t
	
	associate(nu => bandpass(1,1:n), tau => bandpass(2,1:n))
	bandpass_dBdA = integral(tau * nu**(beta+2.0)/sinch(gamma*nu)/exp(gamma*nu), nu, n) * sinch(gamma*mu)*exp(gamma*mu)/mu**(beta+2.0)
	end associate
end function

! modified blackbody color correction, in T[CMB] units
function bbody_log_cc(mu, t, beta, bandpass, n)
	real bbody_log_cc, mu, t, beta, bandpass(:,:); integer n
	
	! calculations are often repeated at the same frequency
	real, save :: dBdT = 0.0, nu0 = -1.0
	
	if (mu /= nu0) dBdT = bandpass_dBdT(mu, 2.7255, 0.0, bandpass, n)
	bbody_log_cc = log(dBdT/bandpass_dBdA(mu, t, beta, bandpass, n)); nu0 = mu
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! magnetic field reconstruction from dust polarization pattern by annealing
subroutine anneal_magnetic(nside, order, map, out, prior)
	real(IO), dimension(0:12*nside**2-1,3) :: map, out
	real(IO), dimension(0:12*nside**2-1), optional :: prior
	integer nside, hside, order, lmax
	
	! allocatable storage for (large) temporary maps
	real(DP), allocatable :: A(:,:,:), B(:,:), L(:,:) ! in NEST ordering !
	real(DP), allocatable :: H(:), W(:), PDF(:)       ! in RING ordering !
	integer, allocatable  :: nn(:,:), idx(:,:)        ! in NEST ordering !
	real(DP), allocatable :: U(:)                     ! temporary vector !
	complex(DPC), allocatable :: alms(:,:,:)
	
	! temporary variables
	real(DP) g(4), v(3), c, f, energy
	integer i, j, k, n, m, ii, step, stride
	
	! magnetic field direction prior [IAU = (70 deg,24 deg), PIP XLIV]
	real(DP), parameter :: B0(3) = [0.31245095, 0.85845193, 0.40673664]
	
	! Allen-Cahn diffusion solver parameters
	real(DP), parameter :: lambda = 2.0/8**2
	real(DP), parameter :: alpha = 0.10, beta = lambda*alpha
	
	! intermediate map filename
	character(len=80) frame
	
	! allocate temporary storage
	hside = max(nside/8, min(nside,8)); lmax = 3*hside-1
	n = nside2npix(nside) - 1; m = nside2npix(hside) - 1
	allocate(A(3,4,0:n), B(0:n,3), U(0:n), H(0:m), W(0:m), PDF(0:m))
	allocate(alms(1,0:lmax,0:lmax), L(9,0:n), nn(9,0:n), idx(4,0:n))
	
	! load pqu map and convert it to NEST oredering
	B = map; if (order == RING) call convert_ring2nest(nside, B)
	
	! initialize possible directions and Laplacian stencil
	do i = 0,n
		A(:,:,i) = pqu2magnetic(nside, NEST, i, B(i,:))
		call stencil(nside, NEST, i, nn(:,i), Lw=L(:,i))
		do k = 1,4; call vec2pix_ring(hside, A(:,k,i), idx(k,i)); end do
	end do
	
	! load prior on magnetic field direction (if any, otherwise use B0)
	if (present(prior)) then
		U = prior; if (order == RING) call convert_ring2nest(nside, U)
		call udgrade_nest(U, nside, W, hside); call convert_nest2ring(nside, W)
	else
		do i = 0,m; call pix2vec_ring(hside, i, v); W(i) = exp(-log(2.0)*(1.0-sum(B0*v))**4); end do
		!call write_map('prior-window.fits', reshape(W,[m+1,1]), hside, RING, creator='FCALC-DEBUG')
	end if
	
	! bootstrap global distribution of magnetic field directions
	PDF = 1.0; do step = 1,4
		H = 0.0; call random_number(U)
		do i = 0,n; k = idx(draw(4, PDF(idx(:,i)), U(i)), i); H(k) = H(k)+1; end do
		
		call map2alm(hside, lmax, lmax, H*W, alms)
		!call alter_alm(hside, lmax, lmax, (20.0*2048)/hside, alms)
		call alm2map(hside, lmax, lmax, alms, PDF)
	end do
	
	! draw a random realization of unit-valued magnetic field
	call random_number(U)
	forall (i=0:n) B(i,:) = A(:,draw(4, PDF(idx(:,i)), U(i)), i)
	
	! minimize energy functional by annealing and diffusion
	do step = 1,640
		! find a coprime stride to traverse the map of (n+1) pixels
		stride = 0; do while (mod(stride,2) == 0 .or. mod(stride,3) == 0)
			call random_number(v); stride = (n+1)/2 * (1.0+v(1)); i = (n+1)*v(2)
		end do
		
		! annealing-diffusion step
		do ii = 0,n; i = i + stride; if (i > n) i = i-(n+1)
			c = norm2(B(i,:))
			
			! flip magnetic field direction to the smoothest one of four possibles
			k = minloc([(sum([(L(j,i)*sum((B(nn(j,i),:)-c*A(:,k,i))**2), j=2,9)]), k=1,4)], 1)
			
			! Allen-Cahn diffusion step on the field projected onto the best direction
			f = alpha*sum(L(2:9,i)*matmul(B(nn(2:9,i),:), A(:,k,i))) + (1.0 + alpha*L(1,i) - beta)*c
			c = f; do j = 1,4; c = (f + (2.0*beta)*c**3)/(1.0 + (3.0*beta)*c**2); end do
			
			! update the field
			B(i,:) = c*A(:,k,i)
		end do
		
		! enforce field average
		B = B/sqrt(sum(B**2)/(n+1))
		
		! energy functional and intermediate maps are output for debugging
		energy = 0.0; do i = 0,n
			energy = energy + sum([(L(j,i)*sum((B(nn(j,i),:)-B(i,:))**2), j=2,9)])/2.0 + (lambda/4.0)*(sum(B(i,:)**2)-1.0)**2
		end do
		
		write (*,*) step, sqrt(sum(B**2)/(n+1)), energy * (pi/3.0)/nside**2
		
		!write (frame,'(g,i0.5,g)') "anneal-", step, ".fits"
		!call write_map(frame, B, nside, NEST, creator='FCALC-DEBUG')
	end do
	
	! copy reconstructed result to output precison/ordering
	out = B; if (order == RING) call convert_nest2ring(nside, out)
	
	deallocate(A, B, U, H, W, PDF, alms, L, nn, idx)
end subroutine

! draw 1 of n values with probabilities p(n)/sum(p)
pure function draw(n, p, v)
	integer i, n, draw; real p(n), c(n), v; intent(in) n, p, v
	
	c(1) = p(1); do i = 2,n; c(i) = c(i-1) + p(i); end do
	do i = 1,n; if (v < c(i)/c(n)) exit; end do; draw = i
end function

end
