#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# whiten kernel data
if [ -r ../cl-lcdm.dat ]; then
    WHITEN=../cl-lcdm.dat
fi

$BINDIR/wiener WHITE:0.0 ${WHITEN:-} > KERNEL-WHITE-00a.dat
$BINDIR/wiener GAUSS:5.0 ${WHITEN:-} > KERNEL-WHITE-05a.dat

for w in 120 800; do
    for k in GAUSS SSG21 SSG42 SSG84; do
       $BINDIR/wiener $k:$w.0 > KERNEL-$k-$w.dat
       $BINDIR/wiener $k:$w.0 ${WHITEN:-} > WHITEN-$k-$w.dat
    done
done
