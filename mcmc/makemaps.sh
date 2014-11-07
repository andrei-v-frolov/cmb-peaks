#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# sims location patern
SIMS=base/cmb-mcmc-????

# create summary directories
if [ ! -d utpmaps ]; then mkdir -pv utpmaps; fi

# remap KS deviation maps (if available)
for f in KS-030 KS-060 KS-090 KS-120 KS-180; do
	MAP="maps/$f.fits"
	
	if [ -f "../$MAP" ]; then
		for p in ltp mtp utp; do
			xTP="utpmaps/$p-${f#KS-}.fits"
			KSD="utpmaps/peaks-$p-ks${f#KS-}"
			NPK="utpmaps/peaks-$p-np${f#KS-}"
			
			echo "Generating log-$p map for $f..."
			$BINDIR/remap $p "../$MAP" "$xTP" $SIMS/$MAP
			
			map2gif -xsz 1200 -bar .true. -sig 2 -inp "$xTP" -out "!$NPK.gif"
			map2gif -xsz 1200 -bar .true. -sig 2 -min 1.3 -max 2.8 -inp "$xTP" -out "!$NPK-95.gif"
			map2gif -xsz 1200 -bar .true. -sig 2 -min 2.0 -max 2.8 -inp "$xTP" -out "!$NPK-99.gif"
		done
		
		map2gif -xsz 1200 -bar .true. -sig 1 -inp "$xTP" -out "!$KSD.gif"
		map2gif -xsz 1200 -bar .true. -sig 1 -min 1.3 -max 2.8 -inp "$xTP" -out "!$KSD-95.gif"
		map2gif -xsz 1200 -bar .true. -sig 1 -min 2.0 -max 2.8 -inp "$xTP" -out "!$KSD-99.gif"
	fi
done
