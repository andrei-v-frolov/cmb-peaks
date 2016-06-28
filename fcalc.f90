! $Id$
! Calculate a pointwise binary operation on two maps
! invoke: fcalc A.fits 'x' B.fits ['=>'] output.fits
! x is a binary operator, see source code for complete list

program fcalc

! HEALPix includes
use mapio
use pdetools
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: lmax = 2000

character(len=8000) :: fin1, op, fin2, fout
integer :: nmaps = 0, nside = 0, ord = 0, n = 0, i
real(IO), dimension(:,:), allocatable :: M1, M2, Mout
logical, dimension(:,:), allocatable :: valid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
call getArgument(1, fin1)
call getArgument(2, op  )
call getArgument(3, fin2)
call getArgument(4, fout)

! syntax sugar
if (fout .eq. '=>') call getArgument(5, fout)

! read input maps
call read_map(fin1, M1, nside, nmaps, ord)
call read_map(fin2, M2, nside, nmaps, ord)
n = nside2npix(nside)-1

! output storage
allocate(Mout, mold=M1); Mout = 0.0

! valid data mask
allocate(valid(0:n,nmaps))
valid = .not. (isnan(M1) .or. isnan(M2))

! apply operator
select case (op)
	! arithmetic operators
	case ('+');  Mout = M1 + M2
	case ('-');  Mout = M1 - M2
	case ('*');  Mout = M1 * M2
	case ('/');  Mout = M1 / M2
	case ('//'); Mout = floor(M1/M2)
	case ('**'); Mout = M1 ** M2
	
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
	case ('inpaint');
		select case (nmaps)
			! case(2) should do tensor inpainting on QU map
			! case(3) should do tensor inpainting on IQU map
			case default; do i = 1,nmaps; call inpaint(M1(:,i), M2(:,i), Mout(:,i), nside, ord); end do
		end select
	
	! polarization operators
	case ('QU->EB');
		select case (nmaps)
			case (2); call rotate_qu2eb(nside, ord, lmax, M1(:,1:2), Mout(:,1:2))
			case (3); call rotate_qu2eb(nside, ord, lmax, M1(:,2:3), Mout(:,2:3)); Mout(:,1) = M1(:,1)
			case default; call abort(trim(op) // " conversion requires QU or IQU map format")
		end select
	
	! unknown operator
	case default; call abort(trim(op) // ": operation not supported")
end select

! write output map
call write_map(fout, Mout, nside, ord, creator='FCALC')

contains

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

end
