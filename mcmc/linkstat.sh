#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# sims location patern
SIMS=base/cmb-mcmc-????

# backup previous results
if [ -f LINKSTAT ]; then mv -f LINKSTAT LINKSTAT.OLD; fi

# create MCMC digests
touch LINKSTAT

for i in $SIMS; do
	echo "python linkstat.py $i/stat/L1-SSG84-*.dat >> LINKSTAT"
done | parallel
