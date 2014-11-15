#!/bin/bash

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=./bins/`uname -m`-`uname -s`

# base kernel to sweep
KERNEL='SSG84'

# sims scaling factor
FUDGE=1.0134

# kernel size sweep
for i in {40..600}; do
	BASE=`printf "L1-$KERNEL-%03i" $i`
	make stat/$BASE.dat gifs/$BASE.gif; rm -f maps/$BASE.fits
	python stat/normalize.py stat/$BASE.dat > stat/$BASE.sig
	#python stat/normalize.py stat/$BASE.dat:mcmc/peakcdf/$BASE.dat:$FUDGE > stat/$BASE.sig
done

# make a merger map
cat stat/L1-$KERNEL-*.sig | $BINDIR/pxl2map maps/white.fits maps/merger.fits
map2gif -inp maps/merger.fits -out '!gifs/merger.gif' -xsz 2400
map2gif -inp maps/merger.fits -out '!gifs/merger-zoom.gif' -pro GNO -lat -56.73 -lon 209.96 -res 3.3 -xsz 1200 -bar true