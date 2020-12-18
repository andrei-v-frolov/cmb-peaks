! Non-linear Wiener filter fitting exponential signal model for foregrounds
! log-wiener lmax {I|IQU}.fits noise-covariance.fits logIQU.fits [residual.fits]

program wiener

! HEALPix includes
use mapio
use extension

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inverse of the (symmetric) noise covariance matrix
pure function inverse(cov)
	real(DP) cov(6), inverse(6), det; intent(in) cov
	
	associate(II => cov(1), IQ => cov(2), IU => cov(3), QQ => cov(4), QU => cov(5), UU => cov(6))
		det = II*QQ*UU - II*QU*QU - IQ*IQ*UU + 2*IQ*IU*QU - IU*IU*QQ
		inverse = [QQ*UU - QU*QU, -IQ*UU + IU*QU, IQ*QU - IU*QQ, II*UU - IU*IU, -II*QU + IQ*IU, II*QQ - IQ*IQ]/det
	end associate
end function

! multiply a 3-vector by a symmetric 3x3 matrix
pure function smatmul(M, a)
	real(DP) M(6), a(3), smatmul(3); intent(in) M, a
	
	associate(II => M(1), IQ => M(2), IU => M(3), QQ => M(4), QU => M(5), UU => M(6), I => a(1), Q => a(2), U => a(3))
		smatmul = [II*I + IQ*Q + IU*U, IQ*I + QQ*Q + QU*U, IU*I + QU*Q + UU*U]
	end associate
end function

! pixel parametrization cost functional and its derivative
pure function cost(iqu, data, invcov)
	real(DP) iqu(3), data(3), invcov(6), cost(4)
	intent(in) iqu, data, invcov
	
	real(DP) p, model(3), y(3)
	
	associate(i => iqu(1), q => iqu(2), u => iqu(3), g => cost(1:3), H => cost(4))
		p = sqrt(q*q + u*u); model = exp(i)*[cosh(p), [q,u]*sinh(p)/p]
		y = smatmul(invcov, data-model); H = dot_product(y, data-model)/2.0
		g = smatmul([model, ([q*q, q*u, u*u]*model(1) + [u*model(3), -u*model(2), q*model(2)])/(p*p)], -y)
	end associate
end function

! non-linear Wiener filter for a single pixel
function fit_iqu(data, invcov, prior)
	real(DP) data(3), invcov(6), prior(3), fit_iqu(3)
	
	! parameters of the optimization problem
	integer, parameter :: n = 3, m = 7, nbd(3) = 0, iprint = -1
	real(DP), parameter :: b(3) = 0.0, eps = 1.0, gtol = 0.0
	
	! optimization variables
	real(DP) f, g(3), x(3), y(4)
	
	! working variables
	character*60 task, csave
	real(DP) wa((2*m+5)*n + 11*m*m + 8*m), dsave(29)
	integer iwa(3*n), isave(44); logical lsave(4)
	
	! initial guess
	task = 'START'; f = norm2(data); x = [log(f), data(2:3)/f]
	
	! optimization loop
	do
		call setulb(n,m,x,b,b,nbd,f,g,eps,gtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)
		if (task(1:5) .eq. 'NEW_X') cycle
		if (task(1:2) .ne. 'FG') exit
		
		y = cost(x, data, invcov)
		g = y(1:3) + prior*x
		f = y(4) + sum(prior*x*x)/2.0
	end do
	
	! check that convergence was achieved
	if (task(1:4) .ne. 'CONV') call abort("optimization failed in fit_iqu(): " // trim(task))
	
	! return best guess
	fit_iqu = x
end function

end
