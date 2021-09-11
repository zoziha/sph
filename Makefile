# Fortran sph Makefile

LIB = sph

FC = gfortran
FFLAGS = -g

export LIB
export FC
export FFLAGS

.PHONY: all clean test

all:
	$(MAKE) -f Makefile --directory=src
	$(MAKE) -f Makefile --directory=app
	$(MAKE) -f Makefile --directory=test

test:
	$(MAKE) -f Makefile --directory=test

clean:
	$(MAKE) -f Makefile clean --directory=src
	$(MAKE) -f Makefile clean --directory=app
	$(MAKE) -f Makefile clean --directory=test