#!/bin/bash

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=./bins/`uname -m`-`uname -s`

# base kernel to sweep
KERNEL='L1-SSG84'

# sims scaling factor
FUDGE=1.0134

# kernel size sweep
for i in {120..999}; do
	make stat/$KERNEL-$i.dat; rm -f maps/$KERNEL-$i.fits
	python stat/normalize.py stat/$KERNEL-$i.dat > stat/$KERNEL-$i.sig
	#python stat/normalize.py stat/$KERNEL-$i.dat:mcmc/peakcdf/$KERNEL-$i.dat:$FUDGE > stat/$KERNEL-$i.sig
done

# make a merger map
cat stat/$KERNEL-*.sig | $BINDIR/pxl2map maps/white.fits maps/merger.fits
map2gif -inp maps/merger.fits -out '!gifs/merger.gif' -xsz 2400
