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
integer i, j, seed(2)

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
	
	! masking operators
	case ('valid'); where (valid) Mout = 1.0
	case ('invalid'); where (.not. valid) Mout = 1.0
	case ('mask'); Mout = M1*M2; where (M2 == 0.0) Mout = 1.0/0.0
	case ('unmask'); Mout = M1/M2; where (M2 == 0.0) Mout = 1.0/0.0
	case ('within:'); where (M1 >= M2 .and. M1 <= M3) Mout = 1.0
	case ('apodize:'); Mout = apodize((M1-M2)/(M3-M2))
	
	! inpainting and filling
	case ('inpaint'); do i = 1,nmaps; call inpaint(nside, ord, M1(:,i), M2(:,i), Mout(:,i)); end do
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
	case ('sum');     nmaps = 1; Mout(:,1) = sum(M1,2)
	case ('product'); nmaps = 1; Mout(:,1) = product(M1,2)
	case ('select');  nmaps = 1; forall (i=0:n) Mout(i,1) = M1(i,M2(i,1))
	
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
			case (1); call randomize(nside, ord, 1, lmax, M1(:,1), Mout(:,1))
			case (3); call randomize(nside, ord, 1, lmax, M1(:,1), Mout(:,1))
			          call randomize_qu(nside, ord, 1, lmax, M1(:,2:3), Mout(:,2:3), .true., .true.)
			case (2); call randomize_qu(nside, ord, 1, lmax, M1(:,1:2), Mout(:,1:2), .true., .true.)
			case default; call abort(trim(op) // " conversion requires I, QU or IQU map format")
		end select
	case ('randomize-blm');
		select case (nmaps)
			case (2); call randomize_qu(nside, ord, 1, lmax, M1(:,1:2), Mout(:,1:2), randomizeB=.true.)
			case (3); call randomize_qu(nside, ord, 1, lmax, M1(:,2:3), Mout(:,2:3), randomizeB=.true.); Mout(:,1) = M1(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	case ('randomize-elm');
		select case (nmaps)
			case (2); call randomize_qu(nside, ord, 1, lmax, M1(:,1:2), Mout(:,1:2), randomizeE=.true.)
			case (3); call randomize_qu(nside, ord, 1, lmax, M1(:,2:3), Mout(:,2:3), randomizeE=.true.); Mout(:,1) = M1(:,1)
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
	
	! reconstruction operators
	case ('magnetic');
		select case (nmaps)
			case (3); call magnetic_fit(nside, ord, 1, 10, M1, Mout)
			case default; call abort(trim(op) // " reconstruction requires IQU/(I+P) map as input")
		end select
	
	! unknown operator
	case default; call abort(trim(op) // ": operation not supported")
end select

! write output map
call write_map(fout, Mout(:,1:nmaps), nside, ord, creator='FCALC')

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
		case ('frac','log','exp','rank','sqrt','valid','invalid','sum','product')
		case ('randomize','shuffle','randomize-alm','randomize-blm','randomize-elm')
		case ('QU->EB','EB->QU')
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
		case ('nest','ring','QU->EB','EB->QU','magnetic')
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
		case ('<','>','<=','>=','=','==','!=','/=','<>')
		case ('project on','orthogonal','accumulate','select')
		case ('valid','invalid','mask','unmask','inpaint','inpaint QU','QU->pure EB','purify')
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
		case ('inpaint'); if (y .eq. 'with' .or. y .eq. 'apodize') ternary = .true.
		case ('accumulate'); if (y .eq. '-') ternary = .true.
		case ('within','apodize'); if (y .eq. ':') ternary = .true.
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
	
	associate (I => iqu(1), Q => iqu(2), U => iqu(3))
		fraction = iqu/(I + sqrt(Q*Q+U*U))
	end associate
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

! exponent of polarization tensor
pure function exp_iqu(iqu)
	real(IO) exp_iqu(3), iqu(3); intent(in) iqu
	real(DP) P, PxP
	
	associate (I => iqu(1), Q => iqu(2), U => iqu(3))
		PxP = Q*Q + U*U; P = sqrt(PxP)
		exp_iqu(1) = exp(I) * cosh(P)
		exp_iqu(2:3) = [Q,U]/P * exp(I) * sinh(P)
	end associate
end function

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

! fit magnetic field potential to polarization fraction map
subroutine magnetic_fit(nside, order, lmin, lmax, map, fit)
	integer nside, order, lmin, lmax
	real(IO), dimension(0:12*nside**2-1,3) :: map, fit
	
	! allocatable storage for (large) temporary maps
	real(DP), allocatable :: pqu(:,:), M1(:,:,:), M2(:,:,:)
	real(DP), allocatable :: field(:,:,:), basis(:,:,:)
	real(DP), allocatable :: A(:,:), B(:,:), pack(:,:)
	integer, allocatable :: pivot(:)
	
	! temporary alms are much smaller, allocate on stack
	real(DP) u, v, chi2(2), lambda
	complex(DP) alms(1,0:lmax,0:lmax)
	integer i, j, k, n, iteration, status
	
	! proposed iterations are tracked in circular buffer
	integer best, next, slow
	
	! allocate temporary storage
	n = nside2npix(nside) - 1; k = (lmax+1)**2 - lmin**2; best = 1; next = 2
	allocate(pqu(0:n,3), M1(0:n,3,2), M2(0:n,3,k), field(0:n,3,2), basis(0:n,3,k))
	allocate(A(k,k), B(k,1), pack(k,2), pivot(k))
	
	! initialize polarization fraction map
	u = minval(map(:,1))
	v = maxval(map(:,1))
	pqu(:,1) = (v-map(:,1))/(v-u)
	pqu(:,2:3) = map(:,2:3)/(v-u)
	if (order == NEST) call convert_nest2ring(nside, pqu)
	
	! initialize magnetic field basis maps
	if (verbose) write (*,*) "Preparing reconstruction basis..."
	
	do i = 1,k; associate(pack => pack(:,1))
		alms = 0.0; pack = 0.0; pack(i) = 1.0
		call unpack_alms(lmin, lmax, pack, alms(1,lmin:lmax,0:lmax))
		call alm2map_magnetic(nside, lmax, lmax, alms, basis(:,:,i))
	end associate; end do
	
	! initial guess for magnetic potential alms
	if (verbose) write (*,*) "Initializing magnetic field guess..."
	
	pack(:,best) = 0.0; chi2(best) = HUGE(chi2); lambda = 0.1; slow = 0
	call map2alm(nside, lmax, lmax, pqu(:,1), alms, [-1.0,1.0])
	call pack_alms(lmin, lmax, alms(1,lmin:lmax,0:lmax), pack(:,next))
	
	if (verbose) write (*,*) "Reconstructing magnetic field potential, RMS(residual):"
	
	do iteration = 1,1000
		! reconstruction residual
		call unpack_alms(lmin, lmax, pack(:,next), alms(1,lmin:lmax,0:lmax))
		call alm2map_magnetic(nside, lmax, lmax, alms, field(:,:,next))
		call magnetic_fit_residual(nside, field(:,:,next), pqu, M1(:,:,next))
		chi2(next) = sum(M1(:,:,next)**2)/(n+1)
		
		if (verbose) write (*,*) sqrt(chi2(next)), sqrt(sum((pack(:,next)-pack(:,best))**2)), lambda
		
		! damping schedule
		if (chi2(next) < chi2(best)) then
			if (chi2(best)/chi2(next) - 1.0 < 1.0e-4) slow = slow + 1
			best = next; if (slow > k) exit
			lambda = lambda/2.0
		else
			lambda = 10.0*lambda + 1.0e-3
		end if
		
		! construct Levenberg–Marquardt matrices
		do i = 1,k
			call magnetic_fit_derivative(nside, field(:,:,best), basis(:,:,i), M2(:,:,i))
			
			B(i,1) = sum(M1(:,:,best)*M2(:,:,i))
			
			do j = 1,i
				A(i,j) = sum(M2(:,:,i)*M2(:,:,j))
				A(j,i) = A(i,j)
			end do
			
			A(i,i) = A(i,i) * (1.0 + lambda)
		end do
		
		! solve for best fit correction
		call dgesv(k, 1, A, k, pivot, B, k, status)
		
		! bail at first sign of trouble
		if (status /= 0) call abort
		
		! proposed update
		next = mod(best,2) + 1
		pack(:,next) = pack(:,best) + B(:,1)
	end do
	
	! cos^2(gamma) map for reconstructed magnetic field
	call unpack_alms(lmin, lmax, pack(:,best), alms(1,lmin:lmax,0:lmax))
	call alm2map_magnetic(nside, lmax, lmax, alms, field(:,:,best))
	forall (i=0:n) fit(i,:) = polarization(field(i,:,best))
	if (order == NEST) call convert_ring2nest(nside, fit)
	
	! clean up allocated storage
	deallocate(pqu, M1, M2, field, basis, A, B, pack, pivot)
end subroutine

! polarization fraction due to magnetic field
pure function polarization(B)
	real(DP) polarization(3), B(3); intent(in) B
	
	polarization = [B(2)**2+B(3)**2,B(3)**2-B(2)**2,-2*B(2)*B(3)]/sum(B**2)
end function

! polarization fraction Jacobian
pure function jacobian(B)
	real(DP) jacobian(3,3), B(3); intent(in) B
	
	jacobian(1,:) = [-(B(2)**2+B(3)**2), B(1)**2, B(1)**2] * B
	jacobian(2,:) = [ (B(2)**2-B(3)**2), -(B(1)**2+2*B(3)**2), (B(1)**2+2*B(2)**2)] * B
	jacobian(3,:) = [ 2*B(1)*B(2)*B(3), -B(3)*(B(1)**2-B(2)**2+B(3)**2), -B(2)*(B(1)**2+B(2)**2-B(3)**2)]
	
	jacobian = jacobian * 2.0/sum(B**2)**2
end function

subroutine magnetic_fit_residual(nside, B, pqu, map)
	integer nside, i, n
	real(DP), dimension(0:12*nside**2-1,3) :: B, pqu, map
	
	n = nside2npix(nside) - 1
	
	forall (i=0:n) map(i,:) = pqu(i,:) - polarization(B(i,:))
end subroutine

subroutine magnetic_fit_derivative(nside, B1, B2, map)
	integer nside, i, n
	real(DP), dimension(0:12*nside**2-1,3) :: B1, B2, map
	
	n = nside2npix(nside) - 1
	
	forall (i=0:n) map(i,:) = matmul(jacobian(B1(i,:)), B2(i,:))
end subroutine

end
