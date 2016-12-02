! power leakage matrix calculated using Wigner 3j symbols
! invoke: leakage mask-cls.fits leakage.fits

program leakage

! HEALPix includes
use imageio
use fwigxjpf

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

integer, parameter :: lmax = 128

real U(0:lmax), W(0:lmax,0:lmax)
real(4) image(1,0:lmax,0:lmax)

integer l1, l2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call fwig_table_init(2*lmax,3)

!$omp parallel
call fwig_temp_init(2*lmax)
!$omp end parallel

U = 1.0

do l1 = 0,lmax; do l2 = 0,lmax
	W(l1,l2) = XXX(l1,l2)
end do; end do

forall (l1=0:lmax,l2=0:lmax) image(1,l1,l2) = log(W(l1,l2)**2/W(l1,l1)/W(l2,l2))/2.0

call image2fits('leakage.fit', image, [0.0,real(lmax)], [0.0,real(lmax)], ['Pi'])

call fwig_temp_free()
call fwig_table_free()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function XXX(l1, l2)
	integer l1, l2, l3; real XXX, W(0:lmax)
	
	!$omp parallel do
	do l3 = 0,lmax
		W(l3) = (2*l3+1)*fwig3jj(2*l1,2*l2,2*l3, 0,0,0)**2
	end do
	
	XXX = (2*l1+1)*(2*l2+1) * sum(W*U)/(4.0*pi)
end function

end program