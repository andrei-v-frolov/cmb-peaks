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
if [ ! -d coldest ]; then mkdir -pv coldest; fi
if [ ! -d hottest ]; then mkdir -pv hottest; fi
if [ ! -d minmaxa ]; then mkdir -pv minmaxa; fi
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
		tail -1 "$i/stat/$k"
	done | sort -k3g > "hottest/$k"
	
	# find min/max average distribution
	echo "Finding min/max avgs in $k"
	if [ -f minmaxa/$k ]; then
		rm -f minmaxa/$k
		touch minmaxa/$k
	fi
	for i in $SIMS; do
		echo "python minmax.py $i/stat/$k >> minmaxa/$k"
	done | parallel
	
	# further analysis for available data
	if [ -f "../stat/$k" ]; then
		# calculate peak distributions
		echo "Calculating peak CDF in $k"
		python resample.py "../stat/$k" $SIMS/"stat/$k" > "peakcdf/$k"
		
		# coldest peak significance
		echo "Calculating outliers in $k"
		echo -n "COLDEST $k " >> SUMMARY
		python bootstrap.py "../stat/$k" "coldest/$k" $FUDGE >> SUMMARY
	fi
done

# plot peak distributions
for k in `head -2 FILES`; do
	echo "Plotting peak CDF in $k..."
	python ../stat/plotpeaks.py "../stat/$k:peakcdf/$k:$FUDGE" ${k%.dat} 10.0
done
