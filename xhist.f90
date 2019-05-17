! $Id$
! HEALPix cross-histogram generator, produces joint distribution of two values in the maps
! xhist M.fits [(min|*):(max|*)/]x-bins[@channel] [(min|*):(max|*)/]y-bins[@channel] output.fits
! xhist M1.fits M2.fits [(min|*):(max|*)/]x-bins[@channel] [(min|*):(max|*)/]y-bins[@channel] output.fits

program xhist

! HEALPix includes
use mapio
use imageio
use extension

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=8000) :: fin1, fin2, fout
integer :: nmaps = 0, nside = 0, ord = 0, pol = -1
real(IO), allocatable :: M1(:,:), M2(:,:), hist(:,:,:)

integer x, y, nx, ny
real(DP) xx(2), yy(2)

integer i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse arguments
select case (nArguments())
	case (4);
		call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord, pol)
		allocate(M2, source=M1)
		call parse_range(2, x, nx, xx, 1, M1)
		call parse_range(3, y, ny, yy, 2, M1)
		call getArgument(4, fout)
	case (5);
		call getArgument(1, fin1); call read_map(fin1, M1, nside, nmaps, ord, pol)
		call getArgument(2, fin2); call read_map(fin2, M2, nside, nmaps, ord, pol)
		call parse_range(3, x, nx, xx, 1, M1)
		call parse_range(4, y, ny, yy, 1, M2)
		call getArgument(5, fout)
	case default; call abort("cannot parse command line arguments")
end select

! initialize accumulators
allocate(hist(3,nx,ny)); hist = 0.0

do i = 0,nside2npix(nside) - 1
	j = floor((M1(i,x)-xx(1))/(xx(2)-xx(1)) * (nx-1)) + 1; if (j < 1 .or. j > nx) cycle
	k = floor((M2(i,y)-yy(1))/(yy(2)-yy(1)) * (ny-1)) + 1; if (k < 1 .or. k > ny) cycle
	
	hist(1,j,k) = hist(1,j,k) + 1
	hist(2,j,k) = hist(2,j,k) + M1(i,x)
	hist(3,j,k) = hist(3,j,k) + M2(i,y)
end do

! write output historgram
call image2fits(fout, hist, xx, yy)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse range specification 
subroutine parse_range(arg, channel, bins, bounds, default, data)
	integer arg, channel, bins, default
	real(DP) bounds(2); real(IO) data(:,:)
	
	integer a, b, c, status
	character(len=80) :: expr
	
	call getArgument(arg, expr)
	
	! parse data channel
	channel = default; c = index(expr,'@',.true.)
	
	if (c > 0) then
		read (expr(c+1:),*,iostat=status) channel
		if (status /= 0) call abort("cannot parse channel in " // trim(expr))
		c = c-1
	else
		c = len(trim(expr))
	end if
	
	! parse number of bins
	bins = 0; b = index(expr(:c),'/')
	
	if (b > 0) then
		read (expr(b+1:c),*,iostat=status) bins
		if (status /= 0) call abort("cannot parse number of bins in " // trim(expr))
		b = b-1
	else
		read (expr(:c),*,iostat=status) bins
		if (status /= 0) call abort("cannot parse number of bins in " // trim(expr))
	end if
	
	! parse histogram range
	a = index(expr(:b),':')
	bounds(1) = minval(data(:,channel))
	bounds(2) = maxval(data(:,channel))
	
	if (a > 0) then
		read (expr(:a-1),*,iostat=status) bounds(1)
		if (status /= 0 .and. expr(:a-1) /= '*') call abort("cannot parse range in " // trim(expr))
		
		read (expr(a+1:b),*,iostat=status) bounds(2)
		if (status /= 0 .and. expr(a+1:b) /= '*') call abort("cannot parse range in " // trim(expr))
	end if
end subroutine

end
