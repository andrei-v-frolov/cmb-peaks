! $Id$
! HEALPix pixel tool; pixels are indexed using existing map format
! 
! output angular separations (in arcmin) between two pixels on a map:
!       pxltool -a[ngdist] <map.fits> <pixel 1> <pixel 2>
! output polar coordinates (in degrees) of a pixel center on a map:
!       pxltool -c[oords] <map.fits> <pixel>
! output pixels on a map containing vector with polar coordinates (in degrees):
!       pxltool -p[ixel] <map.fits> <theta> <phi>

program pxltool

! HEALPix includes
use extension
use head_fits; use healpix_types
use pix_tools; use fitstools

implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=80) :: opt, f, p
integer nmaps, nside, ord, i, j, n
real(DP) :: vi(3), vj(3), t, theta, phi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! input map resolution and ordering
call getArgument(1, opt); call getArgument(2, f)
n = getsize_fits(f, nmaps=nmaps, nside=nside, ordering=ord)

select case (opt(1:2))
        case ('-a') ! angular distance
                call getArgument(3, p); read (p,*) i; call pix2vec(nside, i, vi)
                call getArgument(4, p); read (p,*) j; call pix2vec(nside, j, vj)
                
                call angdist(vi, vj, t)
                
                write (*,*) (180*60/pi) * t
        case ('-c') ! pixel coordinates
                call getArgument(3, p); read (p,*) i; call pix2vec(nside, i, vi)
                
                call vec2ang(vi, theta, phi)
                
                write (*,*) (180/pi) * (/theta, phi/)
        case ('-p') ! pixel index
                call getArgument(3, p); read (p,*) t; theta = (pi/180) * t
                call getArgument(4, p); read (p,*) t;   phi = (pi/180) * t
                
                call ang2vec(theta, phi, vi)
                call vec2pix(nside, vi, i)
                
                write (*,*) i
        case default ! unexpected
                call abort
end select


contains

subroutine pix2vec(nside, p, v)
        integer p, nside; real(DP) v(3)
        
        select case (ord)
                case (1)
                        call pix2vec_ring(nside, p, v)
                case (2)
                        call pix2vec_nest(nside, p, v)
                case default
                        call abort
        end select
end subroutine

subroutine vec2pix(nside, v, p)
        integer p, nside; real(DP) v(3)
        
        select case (ord)
                case (1)
                        call vec2pix_ring(nside, v, p)
                case (2)
                        call vec2pix_nest(nside, v, p)
                case default
                        call abort
        end select
end subroutine
end
