################################################################
# $Id$
# COLDSPOT Pipeline Makefile (Intel Fortran compiler)
################################################################

# Multiple architecture support (sort of)
BINDIR = bins/$(shell uname -m)-$(shell uname -s)
BINS   = $(addprefix $(BINDIR)/,wiener fsynth pxltool pxl2map)

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

BASE  = mask wmap white

KERNS = SSG02 SSG42 SSG84
BEAMS = 060 062 064 065 067 069 071 073 076 078 080 082 085 087 \
	090 092 095 098 101 104 107 110 113 116 120 123 127 131 \
	134 138 142 146 151 155 160 164 169 174 179 184 190 195 \
	201 207 213 219 226 232 239 246 253 260 268 276 284 292 \
	301 309 319 328 337 347 357 368 379 390 401 413 425 437 \
	450 463 477 491 505 520 535 550 566 583 600

FILT = $(shell expr '$*' : '\(\w*\)-[0-9]*')
FWHM = $(shell expr "$*" : '\w*-\([0-9]*\)' \* 2)

MAPS  = $(foreach b,$(BEAMS),$(foreach k,$(KERNS),$(k)-$(b)))

base: $(BINS) $(foreach m,$(BASE),maps/$(m).fits gifs/$(m).gif) alms/mmap.fits alms/mask.fits
maps: base $(foreach m,$(MAPS),maps/$(m)Y.fits)
gifs: maps $(foreach m,$(MAPS),gifs/$(m)Y.gif)
stat: maps $(foreach m,$(MAPS),stat/$(m)Y.dat stat/$(m)Y.pdf)

all: base maps gifs stat

clean:
	rm -f *.bak gmon.out


################### Whiten Pipeline ############################

cl-%.fits alms/%.fits: pars/ana-%.dat maps/wmap.fits
	anafast $<

maps/white.fits: cl-lcdm.dat alms/wmap.fits
	$(BINDIR)/wiener WHITE:0.0 $< | $(BINDIR)/fsynth alms/wmap.fits $@ > /dev/null


################### Convolution Rules ##########################

maps/%B.fits: maps/white.fits
	-lockfile -r0 $@.lock && ($(BINDIR)/$(FILT) $(FWHM).0 $^ $@:3 maps/mask.fits:0.20; rm -f $@.lock)

maps/%Y.fits stat/%Y.dat: alms/mmap.fits alms/mask.fits
	-lockfile -r0 $@.lock && ($(BINDIR)/wiener $(FILT):$(FWHM).0 | $(BINDIR)/fsynth $< maps/$*Y.fits alms/mask.fits:0.5 | sort -k2g > stat/$*Y.dat; rm -f $@.lock)


################### Post-Processing Rules ######################

gifs/%.gif: maps/%.fits
	map2gif -inp $^ -out !$@ -xsz 1200 -bar true

stat/%B.dat: maps/%B.fits
	sed "s/@FILE@/$*/g" < pars/hotspot.dat > /tmp/$*.par
	hotspot /tmp/$*.par; cat /tmp/$*-{min,max}.dat | sort -k2g > $@

stat/%.eps stat/%.pdf: stat/%.dat
	(cd stat; sed "s/@DATA@/$*/g;s/@FILTER@/$(FILT)/g;s/@BEAM@/$(FWHM).0/g" < peak.gpl | gnuplot)


################################################################

$(BINDIR)/wiener: polint.f

$(BINDIR)/%: %.f90
	$(FC) $(FFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)
