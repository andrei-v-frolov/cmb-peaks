################################################################
# $Id$
# COLDSPOT Pipeline Makefile (Intel Fortran compiler)
################################################################

# Multiple architecture support (sort of)
ARCH = $(shell uname -m)-$(shell uname -s)
BINS = bins/$(ARCH)

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

MAPS = mask wmap white
GIFS = $(addprefix gifs/,$(addsuffix .gif,$(MAPS)))

all: $(GIFS)

clean:
	rm -f *.bak gmon.out


################### Whiten Pipeline ############################

cl-wmap.fits alms/wmap.fits: maps/wmap.fits pars/ana-wmap.dat
	anafast pars/ana-wmap.dat

bl-whiten.fits: $(BINS)/whiten cl-lcdm.fits
	$^ $@

alms/white.fits: alms/wmap.fits bl-whiten.fits pars/alt-whiten.dat
	alteralm pars/alt-whiten.dat

maps/white.fits: alms/white.fits pars/syn-whiten.dat
	synfast pars/syn-whiten.dat


################################################################

$(BINS)/%: %.f90
	$(FC) $(FFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)

gifs/%.gif: maps/%.fits
	map2gif -inp $^ -out !$@ -xsz 1200 -bar true
