! $Id$
! Create digest of peak statisitics data (supplied via stdin)

program digest

implicit none

real theta, phi, v, maxs, maxv, mins, minv
integer kind, status, nmax, nmin, n

maxs = 0.0; maxv = 0.0; nmax = 0
mins = 0.0; minv = 0.0; nmin = 0

n = 0

! read data from stdin
do
	read (*,*,iostat=status) theta, phi, v, kind; if (status < 0) exit
	
	if (kind > 0) then; if (v > maxv) maxv = v; maxs = maxs+v; nmax = nmax+1; end if
	if (kind < 0) then; if (v < minv) minv = v; mins = mins+v; nmin = nmin+1; end if
	
	n = n+1
end do

write (*,'(4G24.16,3I)') -mins/nmin, maxs/nmax, -minv, maxv, nmin, nmax, n

end