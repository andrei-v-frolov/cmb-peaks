#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# generate local parameters if not present
if [ ! -f LOCAL-PARAM ]; then
	( cd ..; make stat/L1-GAUSS-020.dat )
	./paramap.py ../stat/L1-GAUSS-020.dat 30.0 > LOCAL-PARAM
fi

# make local parameter maps
awk '{print $1,$2,$3;}' LOCAL-PARAM | $BINDIR/pxl2map order-nest.fits peaks-n.fits
awk '{print $1,$2,$4;}' LOCAL-PARAM | $BINDIR/pxl2map order-nest.fits peaks-ksdev.fits
awk '{print $1,$2,$5;}' LOCAL-PARAM | $BINDIR/pxl2map order-nest.fits peaks-gamma.fits
awk '{print $1,$2,$6;}' LOCAL-PARAM | $BINDIR/pxl2map order-nest.fits peaks-sigma.fits
awk '{print $1,$2,$7;}' LOCAL-PARAM | $BINDIR/pxl2map order-nest.fits peaks-alpha.fits
awk '{print $1,$2,$8;}' LOCAL-PARAM | $BINDIR/pxl2map order-nest.fits peaks-local.fits

# and the corresponding GIFs
map2gif -xsz 800 -bar .true. -inp peaks-n.fits -out '!peaks-n.gif'
map2gif -xsz 800 -bar .true. -inp peaks-ksdev.fits -out '!peaks-ksdev.gif'
map2gif -xsz 800 -bar .true. -min 0.5137003373 -inp peaks-ksdev.fits -out '!peaks-ksdev-95.gif'
map2gif -xsz 800 -bar .true. -min 0.7170542898 -inp peaks-ksdev.fits -out '!peaks-ksdev-99.gif'
map2gif -xsz 800 -bar .true. -inp peaks-gamma.fits -out '!peaks-gamma.gif'
map2gif -xsz 800 -bar .true. -inp peaks-sigma.fits -out '!peaks-sigma.gif'
map2gif -xsz 800 -bar .true. -inp peaks-alpha.fits -out '!peaks-alpha.gif'
map2gif -xsz 800 -bar .true. -inp peaks-local.fits -out '!peaks-local.gif'
map2gif -xsz 800 -bar .true. -min 0.5137003373 -inp peaks-ksdev.fits -out '!peaks-local-95.gif'
map2gif -xsz 800 -bar .true. -min 0.7170542898 -inp peaks-ksdev.fits -out '!peaks-local-99.gif'
