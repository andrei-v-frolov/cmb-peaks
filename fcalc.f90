! $Id$
! HEALPix map calculator, produces output map from zero or more inputs via
!  prefix operator: fcalc 'x' M.fits [=:] output.fits
! postfix operator: fcalc M.fits 'x' [=:] output.fits
!  binary operator: fcalc M1.fits 'x' M2.fits [=:] output.fits
! ternary operator: fcalc M1.fits 'x' M2.fits 'y' M3.fits [=:] output.fits
! see source code for complete list of operators currently implemented

program fcalc

! HEALPix includes
use mapio
use pdetools
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

character(len=8000) :: op, fin1, fin2, fin3, fout
integer :: nmaps = 0, nside = 0, lmax = 0, ord = 0, n = 0
real(IO), dimension(:,:), allocatable :: M1, M2, M3, Mout
logical, dimension(:,:), allocatable :: valid
integer i, seed(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
if (.not. (prefix() .or. postfix() .or. binary() .or. ternary())) call abort("cannot parse command line expression supplied")
if (.not. allocated(M1)) call abort("no map data supplied, I am done")

! map parameters
n = nside2npix(nside)-1; lmax = 3*nside-1

! output storage
allocate(Mout, mold=M1); Mout = 0.0

! valid data mask
allocate(valid(0:n,nmaps), source=.true.)
if (allocated(M1) .and. .not. allocated(M2) .and. .not. allocated(M3)) valid = .not. (isnan(M1))
if (allocated(M1) .and. allocated(M2) .and. .not. allocated(M3)) valid = .not. (isnan(M1) .or. isnan(M2))
if (allocated(M1) .and. allocated(M2) .and. allocated(M3)) valid = .not. (isnan(M1) .or. isnan(M2) .or. isnan(M3))

! initialize random number generator (use urandom on clusters!)
open (333, file="/dev/random", action='read', form='binary')
read (333) seed; call random_seed(PUT=seed); close (333)

! apply operator
select case (op)
	! arithmetic operators
	case ('+');  Mout = M1 + M2
	case ('-');  Mout = M1 - M2
	case ('*');  Mout = M1 * M2
	case ('/');  Mout = M1 / M2
	case ('//'); Mout = floor(M1/M2)
	case ('**'); Mout = M1 ** M2
	case ('accumulate'); Mout = M1 + M2*M2
	case ('accumulate-'); Mout = M1 + (M2-M3)**2
	
	! comparison operators
	case ('<');  where (M1 < M2) Mout = 1.0
	case ('>');  where (M1 > M2) Mout = 1.0
	case ('<='); where (M1 <= M2) Mout = 1.0
	case ('>='); where (M1 >= M2) Mout = 1.0
	case ('=','=='); where (M1 == M2) Mout = 1.0
	case ('!=','/=','<>'); where (M1 /= M2) Mout = 1.0
	
	! projection operators
	case ('project on'); Mout = sum(M1*M2,valid)/sum(M2*M2,valid) * M2
	case ('orthogonal'); Mout = M1 - sum(M1*M2,valid)/sum(M2*M2,valid) * M2
	
	! masking operators
	case ('valid'); where (valid) Mout = 1.0
	case ('invalid'); where (.not. valid) Mout = 1.0
	case ('mask'); Mout = M1*M2; where (M2 == 0.0) Mout = 1.0/0.0
	case ('unmask'); Mout = M1/M2; where (M2 == 0.0) Mout = 1.0/0.0
	
	! inpainting and filling
	case ('inpaint');
		select case (nmaps)
			! case(2) should do tensor inpainting on QU map
			! case(3) should do tensor inpainting on IQU map
			case default; do i = 1,nmaps; call inpaint(M1(:,i), M2(:,i), Mout(:,i), nside, ord); end do
		end select
	case ('inpaintwith');
		select case (nmaps)
			! case(2) should do tensor inpainting on QU map
			! case(3) should do tensor inpainting on IQU map
			case default; do i = 1,nmaps; call inpaint(M1(:,i), M2(:,i), Mout(:,i), nside, ord, fill=M3(:,i)); end do
		end select
	
	! logarithm of a tensor map
	case ('log');
		select case (nmaps)
			case (1); Mout = log(M1)
			case (3); forall (i=0:n) Mout(i,:) = log_iqu(M1(i,:))
			case default; call abort(trim(op) // " conversion requires I or IQU map format")
		end select
	
	! random map generators
	case ('randomize');
		allocate(M2, mold=M1); allocate(M3, mold=M1)
		call random_number(M2); call random_number(M3)
		Mout = M1 * sqrt(-2.0*log(M2)) * cos(2.0*pi*M3)
	
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
	
	! unknown operator
	case default; call abort(trim(op) // ": operation not supported")
end select

! write output map
call write_map(fout, Mout, nside, ord, creator='FCALC')

contains

! parse prefix operator command line
function prefix()
	character(len=80) :: x; logical prefix; prefix = .false.
	
	! argument count guard
	if (nArguments() < 3 .or. nArguments() > 4) return
	
	! operator placement
	call getArgument(1, x)
	
	! prefix operation guard
	select case (x)
		case ('log','valid','invalid','randomize','QU->EB','EB->QU')
		case default; return
	end select
	
	! operator name
	prefix = .true.; op = trim(x)
	
	! read input maps
	call getArgument(2, fin1); call read_map(fin1, M1, nside, nmaps, ord)
	
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
		case ('QU->EB','EB->QU')
		case default; return
	end select
	
	! operator name
	postfix = .true.; op = trim(x)
	
	! read input maps
	call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord)
	
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
		case ('project on','orthogonal', 'accumulate')
		case ('<','>','<=','>=','=','==','!=','/=','<>')
		case ('valid','invalid','mask','unmask','inpaint')
		case default; return
	end select
	
	! operator name
	binary = .true.; op = trim(x)
	
	! read input maps
	call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord)
	call getArgument(3, fin2); call read_map(fin2, M2, nside, nmaps, ord)
	
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
		case ('inpaint'); if (y .eq. 'with') ternary = .true.
		case ('accumulate'); if (y .eq. '-') ternary = .true.
	end select
	
	! ternary operation guard
	if (.not. ternary) return; op = trim(x) // trim(y)
	
	! read input maps
	call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord)
	call getArgument(3, fin2); call read_map(fin2, M2, nside, nmaps, ord)
	call getArgument(5, fin3); call read_map(fin3, M3, nside, nmaps, ord)
	
	! output map name
	call getArgument(6, fout); if (fout .eq. '=:') call getArgument(7, fout)
end function

! logarithm of polarization tensor
pure function log_iqu(iqu)
	real(IO) log_iqu(3), iqu(3); intent(in) iqu
	real(DP) P, PxP
	
	associate (I => iqu(1), Q => iqu(2), U => iqu(3))
		PxP = Q*Q + U*U; P = sqrt(PxP)
		log_iqu(1) = log(I*I - PxP)/2.0
		log_iqu(2:3) = [Q,U]/P * log((I+P)/(I-P))/2.0
	end associate
end function

! full-sky QU to EB rotation wrapper
subroutine rotate_qu2eb(nside, order, lmax, QU, EB)
	use alm_tools
	
	integer nside, npix, lmax, order, spin
	real(IO), dimension(0:12*nside**2-1,1:2) :: QU, EB
	real(DP), dimension(:,:), allocatable :: map
	complex(DPC), allocatable :: alms(:,:,:)
	
	npix = nside2npix(nside); spin = 2
	allocate(map(0:npix-1,1:2), alms(1:2, 0:lmax, 0:lmax))
	
	map = QU; if (order == NEST) call convert_nest2ring(nside, map)
	
	call map2alm_spin(nside, lmax, lmax, spin, map, alms)
	call alm2map(nside, lmax, lmax, alms(1:1,:,:), map(:,1))
	call alm2map(nside, lmax, lmax, alms(2:2,:,:), map(:,2))
	
	if (order == NEST) call convert_ring2nest(nside, map); EB = map
	
	deallocate(map, alms)
end subroutine

! full-sky EB to QU rotation wrapper
subroutine rotate_eb2qu(nside, order, lmax, EB, QU)
	use alm_tools
	
	integer nside, npix, lmax, order, spin
	real(IO), dimension(0:12*nside**2-1,1:2) :: EB, QU
	real(DP), dimension(:,:), allocatable :: map
	complex(DPC), allocatable :: alms(:,:,:)
	
	npix = nside2npix(nside); spin = 2
	allocate(map(0:npix-1,1:2), alms(1:2, 0:lmax, 0:lmax))
	
	map = EB; if (order == NEST) call convert_nest2ring(nside, map)
	
	call map2alm(nside, lmax, lmax, map(:,1), alms(1:1,:,:))
	call map2alm(nside, lmax, lmax, map(:,2), alms(2:2,:,:))
	call alm2map_spin(nside, lmax, lmax, spin, alms, map)
	
	if (order == NEST) call convert_ring2nest(nside, map); QU = map
	
	deallocate(map, alms)
end subroutine

end
