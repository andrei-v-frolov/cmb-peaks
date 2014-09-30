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

character(len=80) :: fin1, op, fin2, fout
integer :: nmaps = 0, nside = 0, ord = 0, n = 0
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
	case ('inpaint'); call inpaint(M1, M2, Mout, nside, ord)
	
	! unknown operator
	case default; call abort(trim(op) // ": operation not supported")
end select

! write output map
call write_map(fout, Mout, nside, ord, creator='FCALC')

end
