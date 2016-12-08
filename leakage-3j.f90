! power leakage matrix calculated using Wigner 3j symbols
! invoke: leakage-3j mask-cls.fits leakage.fits

program leakage_3j

! HEALPix includes
use imageio
use fwigxjpf
use fitstools
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

character(len=8000) :: mask, output
character(len=80) :: header(64), units(64)
real(DP), allocatable :: Cl(:,:), W(:,:,:)

integer, parameter :: lcut = 128
integer lmax, ncl, np, l1, l2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
if (nArguments() /= 2) call abort("cannot parse command line expression supplied")
call getArgument(1, mask); call getArgument(2, output)

! spectrum parameters
np = getsize_fits(mask, nmaps=ncl, mlpol=lmax); if (lmax > lcut) lmax = lcut
if (ncl > 1) write (0,*) "warning: only temperature mask Cl's will be used in leakage calculation!"

! allocate storage and read in mask spectra
allocate(Cl(0:lmax, 1:ncl), W(4,0:lmax,0:lmax))
call fits2cl(mask, cl, lmax, ncl, header, units)

! initialze Wigner 3j tables
call fwig_table_init(2*lmax,3)

!$omp parallel
call fwig_temp_init(2*lmax)
!$omp end parallel

! compute leakage matrix
!$omp parallel do
do l1 = 0,lmax; do l2 = 0,lmax
	W(:,l1,l2) = leak(l1,l2)
end do; end do

!forall (l1=0:lmax,l2=0:lmax) W(1,l1,l2) = log(W(1,l1,l2)**2/W(1,l1,l1)/W(1,l2,l2))/2.0

call image2fits(output, W, [0.0,real(lmax)], [0.0,real(lmax)], ['K','K+','K-','Kx'])

call fwig_temp_free()
call fwig_table_free()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function leak(l1, l2)
	integer l1, l2, l3; real J0, J2, leak(4)
	
        ! initialize accumulators
        leak = 0.0
        
        ! even l1+l2+l3 terms
	do l3 = mod(l1+l2,2),lmax,2
                J0 = fwig3jj(2*l1,2*l2,2*l3,   0,    0,   0)
                J2 = fwig3jj(2*l1,2*l2,2*l3, 2*2,2*(-2),  0)
                leak = leak + (2*l3+1)*Cl(l3,1) * [J0*J0, J2*J2, 0.0, J0*J2]
	end do

        ! odd l1+l2+l3 terms
        do l3 = mod(l1+l2+1,2),lmax,2
                J0 = fwig3jj(2*l1,2*l2,2*l3,   0,    0,   0)
                J2 = fwig3jj(2*l1,2*l2,2*l3, 2*2,2*(-2),  0)
                leak = leak + (2*l3+1)*Cl(l3,1) * [J0*J0, 0.0, J2*J2, 0.0]
        end do
	
	leak = (2*l1+1)*(2*l2+1)/(4.0*pi) * leak/[1,2,2,1]
end function

end program
