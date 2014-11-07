#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# sims location patern
SIMS=base/cmb-mcmc-????

# remap KS deviation maps (if available)
for f in KS-030 KS-060 KS-090 KS-120 KS-180; do
	MAP="maps/$f.fits"
	
	if [ -f "../$MAP" ]; then
		for p in ltp mtp utp; do
			xTP="$p-${f#KS-}.fits"
			KSD="peaks-$p-ks${f#KS-}"
			NPK="peaks-$p-np${f#KS-}"
			
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

# create summary directories
if [ ! -d coldest ]; then mkdir -pv coldest; fi
if [ ! -d hottest ]; then mkdir -pv hottest; fi
if [ ! -d peakcdf ]; then mkdir -pv peakcdf; fi

# backup previous results
if [ -f SUMMARY ]; then mv -f SUMMARY SUMMARY.OLD; fi

# create MCMC digests
for k in `cat FILES`; do
	# find coldest peak distribution
	echo "Finding coldest peak in $k"
	for i in $SIMS; do
		head -1 "$i/stat/$k"
	done | sort -k3g > "coldest/$k"
	
	# find hottest peak distribution
	echo "Finding hottest peak in $k"
	for i in $SIMS; do
		head -1 "$i/stat/$k"
	done | sort -k3g > "hottest/$k"
	
	# further analysis for available data
	if [ -f "../stat/$k" ]; then
		# calculate peak distributions
		echo "Calculating peak CDF in $k"
		python resample.py "../stat/$k" $SIMS/"stat/$k" > "peakcdf/$k"
		
		# coldest peak significance
		echo "Calculating outliers in $k"
		echo -n "COLDEST $k " >> SUMMARY
		python bootstrap.py "../stat/$k" "coldest/$k" >> SUMMARY
	fi
done
