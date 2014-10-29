#!/bin/bash

# abort on errors
set -u -e -o errexit

# binaries location
BINS=$BASE/bins/`uname -m`-`uname -s`

make
make stat/L1-GAUSS-020.dat
make stat/L1-SSG84-400.dat
make maps
make clean

# make plots
(cd kern; ./kernels.sh || true)
(cd xtra; ./paramap.sh || true)
(cd sign; ./plotsign.sh || true)
(cd mcmc; ./summarize.sh || true)

