#!/bin/sh

# abort on errors
set -u -e -o errexit

# binaries location
BINDIR=../bins/`uname -m`-`uname -s`

# sims location patern
SIMS=base/cmb-mcmc-????

# sims scaling factor
FUDGE=1.0134

# create summary directories
if [ ! -d ksdistr ]; then mkdir -pv ksdistr; fi

# backup previous results
if [ -f SUMMARY ]; then mv -f SUMMARY SUMMARY.OLD; fi

# create MCMC digests
for k in `cat FILES`; do
	# find Kolmogorov-Smirnov distribution
	if [ -f "peakcdf/$k" ]; then
		echo "Finding Kolmogorov-Smirnov statistics for $k"
		for i in $SIMS; do
			echo "python ks-test.py $i/stat/$k peakcdf/$k"
		done | parallel | sort -k1g > "ksdistr/$k"
	fi
	
	# calculate Kolmogorov-Smirnov significance
	if [ -f "../stat/$k" ]; then
		echo "Calculating Kolmogorov-Smirnov significance for $k"
		echo -n "KSDEVSG $k " >> SUMMARY
		python ks-test.py "../stat/$k" "peakcdf/$k" $FUDGE | python bootstrap.py logutp "-:0" "ksdistr/$k:0" >> SUMMARY
	fi
done
