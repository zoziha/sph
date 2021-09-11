.SUFFIXES: .f90 .o

SRCF = \
	art_heat.f \
	art_visc.f \
	av_vel.f \
	density.f \
	direct_find.f \
	eos.f \
	external_force.f \
	grid_geom.f \
	hsml.f \
	init_grid.f \
	input.f \
	internal_force.f \
	kernel.f \
	link_list.f \
	output.f \
	single_step.f \
	time_integration.f \
	virt_part.f \
	viscosity.f

SRCF90 = \
	time_elapsed.f90 \
	time_print.f90

OBJF   := $(SRCF:.f=.o)
OBJF90 := $(SRCF90:.f90=.o)

lib$(LIB).a: $(OBJF) $(OBJF90)
	ar -rcs lib$(LIB).a $(OBJF) $(OBJF90)

clean:
	rm -f -r *.o *.a *.so *.mod *.smod

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

time_elapsed.o: time_elapsed.f90
time_print.o  : time_print.f90