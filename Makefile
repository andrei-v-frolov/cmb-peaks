################################################################
# $Id$
# COLDSPOT Pipeline Makefile (Intel Fortran compiler)
################################################################

# Multiple architecture support (sort of)
BINDIR = bins/$(shell uname -m)-$(shell uname -s)
BINS   = $(addprefix $(BINDIR)/,wiener fsynth lmask fdiff pxltool pxl2map)

# Fortran compiler (adjust for your machine, -r8 is mandatory)
FC = ifort
FFLAGS = -O3 -ipo -xT -r8 -pc80 -parallel
LDFLAGS = -static-intel

# HEALPix and CFITSIO libraries
HPXINC = -I$(HEALPIX)/include
HPXLIB = -L$(HEALPIX)/lib -lhealpix
FITLIB = -L/opt/local/lib -lcfitsio

INCS += $(HPXINC) $(FITINC)
LIBS += $(HPXLIB) $(FITLIB)


################################################################

BASE := wmap white mask $(foreach k,1 2 3 4,lmap-$(k))
BASE := $(foreach m,$(BASE),maps/$(m).fits alms/$(m).fits gifs/$(m).gif)
BASE := $(filter-out maps/wmap.fits maps/mask.fits alms/white.fits,$(BASE))

STATS = L1 L2 L3 L4 T3 T4
KERNS = SSG02 SSG42 SSG84
BEAMS = 060 062 064 065 067 069 071 073 076 078 080 082 085 087 \
	090 092 095 098 101 104 107 110 113 116 120 123 127 131 \
	134 138 142 146 151 155 160 164 169 174 179 184 190 195 \
	201 207 213 219 226 232 239 246 253 260 268 276 284 292 \
	301 309 319 328 337 347 357 368 379 390 401 413 425 437 \
	450 463 477 491 505 520 535 550 566 583 600

# portable regexp definitions
DIGITS = [0-9]\{1,\}
ALNUMS = [[:alnum:]\-]\{1,\}

FILT = $(shell expr '$*' : '\($(ALNUMS)\)-$(DIGITS)')
FWHM = $(shell expr '$*' : '$(ALNUMS)-\($(DIGITS)\)' \* 2)
KERN = $(shell expr '$*' : '$(ALNUMS)-\($(ALNUMS)\)-$(DIGITS)')

MAPS  = $(foreach b,$(BEAMS),$(foreach k,$(KERNS),$(foreach s,$(STATS),$(s)-$(k)-$(b))))

base: $(BINS) $(BASE)
maps: base $(foreach m,$(MAPS),maps/$(m).fits)
gifs: base $(foreach m,$(MAPS),gifs/$(m).gif)
stat: base $(foreach m,$(filter L1-%,$(MAPS)),stat/$(m).dat stat/$(m).pdf)

all: base maps gifs stat

clean:
	rm -f $(foreach m,$(MAPS),maps/$(m).fits)

dataclean:
	rm -f $(foreach m,$(MAPS),maps/$(m).fits gifs/$(m).gif stat/$(m).{dat,eps,pdf})

distclean: dataclean
	rm -f $(BINS) $(BASE) cl-{wmap,mask}.fits cl-lmap-?.fits


################### Whiten Pipeline ############################

pars/ana-lmap-%.dat: pars/ana-lmap.dat
	sed "s/@ORDER@/$*/g" < $< > $@

cl-%.fits alms/%.fits: pars/ana-%.dat maps/%.fits
	anafast $<

maps/white.fits: cl-lcdm.dat alms/wmap.fits
	$(BINDIR)/wiener WHITE:0.0 $< | $(BINDIR)/fsynth alms/wmap.fits $@ > /dev/null

maps/lmask.fits $(foreach k, 1 2 3 4, maps/lmap-$(k).fits): maps/white.fits maps/mask.fits
	$(BINDIR)/lmask $< $@:4 maps/mask.fits maps/lmap


################### Convolution Rules ##########################

lock = -lockfile -r0 $@.lock && ($(1); rm -f $@.lock)
convolve = $(BINDIR)/fsynth $< $(2) alms/mask.fits:0.5 < $(1) | sort -k2g > $(3)

beam/%.dat:
	$(call lock,$(BINDIR)/wiener $(FILT):$(FWHM).0 > $@)

maps/L1-%.fits stat/L1-%.dat: alms/lmap-1.fits alms/mask.fits beam/%.dat
	$(call lock,$(call convolve,beam/$*.dat,maps/L1-$*.fits,stat/L1-$*.dat))

maps/L2-%.fits stat/L2-%.dat: alms/lmap-2.fits alms/mask.fits beam/%.dat
	$(call lock,$(call convolve,beam/$*.dat,maps/L2-$*.fits,stat/L2-$*.dat))

maps/L3-%.fits stat/L3-%.dat: alms/lmap-3.fits alms/mask.fits beam/%.dat
	$(call lock,$(call convolve,beam/$*.dat,maps/L3-$*.fits,stat/L3-$*.dat))

maps/L4-%.fits stat/L4-%.dat: alms/lmap-4.fits alms/mask.fits beam/%.dat
	$(call lock,$(call convolve,beam/$*.dat,maps/L4-$*.fits,stat/L4-$*.dat))

maps/T3-%.fits: maps/L3-%.fits maps/L2-%.fits
	$(call lock,$(BINDIR)/fdiff $^ $@)

maps/T4-%.fits: maps/L4-%.fits maps/L2-%.fits
	$(call lock,$(BINDIR)/fdiff $^ $@)


################### Post-Processing Rules ######################

gifs/%.gif: maps/%.fits
	map2gif -inp $^ -out !$@ -xsz 1200 -bar true

stat/%.eps stat/%.pdf: stat/%.dat
	(cd stat; sed "s/@DATA@/$*/g;s/@KERNEL@/$(KERN)/g;s/@FILTER@/$(FILT)/g;s/@BEAM@/$(FWHM).0/g" < peak.gpl | gnuplot)


################### Binaries & Dependencies ####################

$(BINDIR)/lmask: rank.f
$(BINDIR)/wiener: polint.f

$(BINDIR)/%: %.f90
	$(FC) $(FFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)


################### Brute Force (obsolete) #####################

maps/%B.fits: maps/white.fits
	$(call lock,$(BINDIR)/$(FILT) $(FWHM).0 $^ $@:3 maps/mask.fits:0.20)

stat/%B.dat: maps/%B.fits
	sed "s/@FILE@/$*/g" < pars/hotspot.dat > /tmp/$*.par
	hotspot /tmp/$*.par; cat /tmp/$*-{min,max}.dat | sort -k2g > $@
