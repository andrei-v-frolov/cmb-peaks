################################################################
# COLDSPOT Pipeline Makefile (Intel Fortran compiler)
################################################################

# Multiple architecture support (sort of)
BINDIR := bins/$(shell uname -m)-$(shell uname -s)
BINS   := $(addprefix $(BINDIR)/,wiener fsynth lmask fcalc leakage-3j leakage-mc pxl2map digest ksmap remap)

# Fortran compiler (adjust for your machine, -r8 is mandatory)
FC = ifort
FFLAGS = -O3 -ipo -xHOST -fpp -heap-arrays 256 -r8 -pc80 -parallel
FFLAGS_OMP = $(subst parallel,qopenmp,$(FFLAGS))
#LDFLAGS = -static-intel

# HEALPix libraries
HEALPIX ?= /opt/healpix
HPXINC = -I$(HEALPIX)/include
HPXLIB = -L$(HEALPIX)/lib -lhealpix

# CFITSIO libraries
CFITSIO_DIR ?= /opt/local
FITINC = -I$(CFITSIO_DIR)/include
FITLIB = -L$(HEALPIX)/lib -lcfitsio -lcurl

# WIGXJPF libraries
WIGXJPF ?= /opt/wigxjpf
WIGINC = -I$(WIGXJPF)/mod
WIGLIB = -L$(WIGXJPF)/lib -lwigxjpf

# LAPACK libraries (use MKL if compiling with Intel Fortran)
MKLROOT ?= /opt/intel/mkl
LAPINC = -I$(MKLROOT)/include
LAPLIB = -L$(MKLROOT)/lib -lmkl_rt

# Intel's dynamic libraries location
LDFLAGS += -Wl,-rpath,/opt/intel/lib
LDFLAGS += -Wl,-rpath,$(MKLROOT)/lib

INCS += $(HPXINC) $(FITINC) $(LAPINC)
LIBS += $(HPXLIB) $(FITLIB) $(LAPLIB)


################################################################

BASE := wmap white mask $(foreach k,1 2 3 4,lmap-$(k))
BASE := $(foreach m,$(BASE),maps/$(m).fits alms/$(m).fits gifs/$(m).gif)
BASE := $(filter-out maps/wmap.fits maps/mask.fits alms/white.fits,$(BASE))

#STATS = L1 L2 L3 L4 T3 T4
#KERNS = SSG02 SSG42 SSG84

STATS = L1
KERNS = GAUSS SSG21 SSG42 SSG84
BEAMS = 060 062 064 065 067 069 071 073 076 078 080 082 085 087 \
	090 092 095 098 101 104 107 110 113 116 120 123 127 131 \
	134 138 142 146 151 155 160 164 169 174 179 184 190 195 \
	201 207 213 219 226 232 239 246 253 260 268 276 284 292 \
	301 309 319 328 337 347 357 368 379 390 401 413 425 437 \
	450 463 477 491 505 520 535 550 566 583 600

# portable regexp definitions
DIGITS = [0-9]\{1,\}
ALNUMS = [[:alnum:]\-]\{1,\}

EXPR := $(word 1,$(shell which gexpr expr))
FILT = $(shell $(EXPR) '$*' : '\($(ALNUMS)\)-$(DIGITS)')
FWHM = $(shell $(EXPR) '$*' : '$(ALNUMS)-\($(DIGITS)\)' \* 2)
KERN = $(shell $(EXPR) '$*' : '$(ALNUMS)-\($(ALNUMS)\)-$(DIGITS)')

MAPS  = $(foreach b,$(BEAMS),$(foreach k,$(KERNS),$(foreach s,$(STATS),$(s)-$(k)-$(b))))

base: $(BINS) $(BASE)
beam: base $(foreach b,$(BEAMS),$(foreach k,$(KERNS),beam/$(k)-$(b).dat))
maps: base $(foreach m,$(MAPS),maps/$(m).fits)
gifs: base $(foreach m,$(MAPS),gifs/$(m).gif)
stat: base $(foreach m,$(filter L1-%,$(MAPS)),stat/$(m).dat stat/$(m).pdf)

all: base beam maps gifs stat

clean:
	rm -f $(foreach m,$(MAPS),maps/$(m).fits)

dataclean:
	rm -f $(foreach m,$(MAPS),maps/$(m).fits gifs/$(m).gif stat/$(m).{dat,eps,pdf})
	rm -f kern/*.{dat,pdf} xtra/peaks-*.{fits,gif}

distclean: dataclean
	rm -f $(BINS) $(BASE) cl-{wmap,mask}.fits cl-lmap-?.fits
	rm -f `find . -name "*.o" -or -name "*.mod" -or -name "*.py[cod]"`


################### Whiten Pipeline ############################

pars/ana-lmap-%.dat: pars/ana-lmap.dat
	sed "s/@ORDER@/$*/g" < $< > $@

cl-%.fits alms/%.fits: pars/ana-%.dat maps/%.fits
	anafast $<

maps/white.fits: cl-lcdm.dat alms/wmap.fits
	$(BINDIR)/wiener GAUSS:5.0 $< | $(BINDIR)/fsynth alms/wmap.fits $@ > /dev/null

maps/lmask.fits $(foreach k, 1 2 3 4, maps/lmap-$(k).fits): maps/white.fits maps/mask.fits
	$(BINDIR)/lmask $< $@:4 maps/mask.fits maps/lmap


################### Convolution Rules ##########################

lock = -ln -s $@ $@.lock && ($(1); rm -f $@.lock)
convolve = $(BINDIR)/fsynth $< $(2) alms/mask.fits:0.5 < $(1) | sort -k3g > $(3)

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
	$(call lock,$(BINDIR)/fcalc $(word 1,$^) '/' $(word 2,$^) '=>' $@)

maps/T4-%.fits: maps/L4-%.fits maps/L2-%.fits
	$(call lock,$(BINDIR)/fcalc $(word 1,$^) '/' $(word 2,$^) '=>' $@)


################### Post-Processing Rules ######################

gifs/%.gif: maps/%.fits
	map2gif -inp $^ -out !$@ -xsz 1200 -bar true
#	map2gif -inp $^ -out !$@ -pro GNO -lat -56.732571412162 -lon 209.958217270195 -res 3.3 -xsz 800 -bar true

stat/%.pdf: stat/%.dat
	-(cd stat; python plotpeaks.py $*.dat $* 12.0)


################### Binaries & Dependencies ####################

# binaries
$(BINDIR)/fcalc: mapio.o pdetools.o almtools.o rank.o
$(BINDIR)/fsynth: mapio.o rank.o
$(BINDIR)/lmask: mapio.o rank.o
$(BINDIR)/ksmap: mapio.o rank.o
$(BINDIR)/remap: mapio.o rank.o
$(BINDIR)/pxl2map: mapio.o
$(BINDIR)/wiener: polint.o

# OpenMP binaries
$(BINDIR)/leakage-mc: mapio.o imageio.o pdetools.o almtools.o
$(BINDIR)/leakage-3j: leakage-3j.f90 imageio.o
	$(FC) $(FFLAGS_OMP) $(WIGINC) $(INCS) $^ -o $@ $(LDFLAGS) $(WIGLIB) $(LIBS)

# modules
mapio.o: mapio.fin
imageio.o: imageio.fin
pdetools.o: multigrid.fin inpaint-qu.fin almtools.o
almtools.o: almtools.fin maptools.fin complex-qu.fin mapio.o

# generic rules
$(BINDIR)/%: %.f90
	$(FC) $(FFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)
%.o %.mod: %.f90
	$(FC) $(FFLAGS) $(INCS) $< -c
