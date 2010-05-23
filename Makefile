################################################################
# $Id$
# COLDSPOT Pipeline Makefile (Intel Fortran compiler)
################################################################

# Fortran compiler (adjust for your machine, -r8 is mandatory)
FC = ifort
FFLAGS = -O3 -ipo -xT -r8 -pc80
LDFLAGS = -static-intel

# HEALPix library
HPXINC = -I$(HEALPIX)/include
HPXLIB = -L$(HEALPIX)/lib -lhealpix -lcfitsio


################################################################

MAPS = mask wmap white
GIFS = $(addprefix gifs/,$(addsuffix .gif,$(MAPS)))

all: $(GIFS)

clean:
	rm -f *.bak gmon.out

cl-wmap.fits alms/wmap.fits: maps/wmap.fits pars/ana-wmap.dat
	anafast pars/ana-wmap.dat

bl-whiten.fits: whiten cl-lcdm.fits
	$^ $@

alms/white.fits: alms/wmap.fits bl-whiten.fits pars/alt-whiten.dat
	alteralm pars/alt-whiten.dat

maps/white.fits: alms/white.fits pars/syn-whiten.dat
	synfast pars/syn-whiten.dat


################################################################

%: %.f90
	$(FC) $(FFLAGS) $(HPXINC) $^ -o $@ $(LDFLAGS) $(HPXLIB)

gifs/%.gif: maps/%.fits
	map2gif -inp $^ -out !$@ -xsz 1200 -bar true
