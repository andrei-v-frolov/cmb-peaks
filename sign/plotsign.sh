#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# make coldspot significance files
for k in GAUSS SSG21 SSG42 SSG84; do
	python coldspots.py ../stat/L1-$k-*.dat > $k.dat
done
