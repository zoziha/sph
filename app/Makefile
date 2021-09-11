all: sph

# Orginal test
sph: sph.o 
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	cd .. && time ./app/sph.x