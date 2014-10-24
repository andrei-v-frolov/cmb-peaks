#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# sims location patern
SIMS=base/cmb-mcmc-????

# remap local deviation maps if available
for f in peaks-ksdev peaks-local; do
	MAP="xtra/$f.fits"
	UTP="peaks-utp-${f#peaks-}.fits"
	
	if [ -f "../$MAP" ]; then
		echo "Generating log-UTP map for $f..."
		$BINDIR/remap "../$MAP" "$UTP" $SIMS/$MAP
		map2gif -xsz 800 -bar .true. -inp "$UTP" -out "!${UTP%.fits}.gif"
		map2gif -xsz 800 -bar .true. -min 1.3 -max 2.7 -inp "$UTP" -out "!${UTP%.fits}-95.gif"
		map2gif -xsz 800 -bar .true. -min 2.0 -max 2.7 -inp "$UTP" -out "!${UTP%.fits}-99.gif"
	fi
done

# create summary directories
if [ ! -d coldest ]; then mkdir -pv coldest; fi
if [ ! -d hottest ]; then mkdir -pv hottest; fi

# backup previous results
if [ -f SUMMARY ]; then mv -f SUMMARY SUMMARY.OLD; fi

# find coldest peak distribution
for k in `cat FILES`; do
	echo "Finding coldest peak in $k"
	for i in $SIMS; do
		head -1 "$i/stat/$k"
	done | sort -k3g > "coldest/$k"
	if [ -f "../stat/$k" ]; then
		echo -n "COLDEST $k " >> SUMMARY
		python bootstrap.py "../stat/$k" "coldest/$k" >> SUMMARY
	fi
done

# find hottest peak distribution
for k in `cat FILES`; do
	echo "Finding hottest peak in $k"
	for i in $SIMS; do
		head -1 "$i/stat/$k"
	done | sort -k3g > "hottest/$k"
done
