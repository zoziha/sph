
!---------------------------------------------------------
!     including file for parameters and constants used
!     in the entire sph software packages.
!---------------------------------------------------------

!     dim : dimension of the problem (1, 2 or 3)
    integer :: dim
    parameter(dim=2)

!     maxn    : maximum number of particles (gfortran: *<WARNING>* `-fmax-stack-var-size=`, r change the code to use an ALLOCATABLE array)
!     max_interation : maximum number of interaction pairs
    integer :: maxn, max_interaction
    parameter(maxn=12000, max_interaction=100*maxn)

!     parameters for the computational geometry,
!     x_maxgeom : upper limit of allowed x-regime
!     x_mingeom : lower limit of allowed x-regime
!     y_maxgeom : upper limit of allowed y-regime
!     y_mingeom : lower limit of allowed y-regime
!     z_maxgeom : upper limit of allowed z-regime
!     z_mingeom : lower limit of allowed z-regime
    real(rk) :: x_maxgeom, x_mingeom, y_maxgeom, y_mingeom, z_maxgeom, z_mingeom
    parameter(x_maxgeom=10._rk, x_mingeom=-10._rk, y_maxgeom=10._rk, y_mingeom=-10._rk, z_maxgeom=10._rk, z_mingeom=-10._rk)

!     sph algorithm for particle approximation (pa_sph)
!     pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
!              2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
    integer :: pa_sph
    parameter(pa_sph=2)

!     nearest neighbor particle searching (nnps) method
!     nnps = 1 : simplest and direct searching
!            2 : sorting grid linked list
!            3 : tree algorithm
    integer :: nnps
    parameter(nnps=1)

!     smoothing length evolution (sle) algorithm
!     sle = 0 : keep unchanged,
!           1 : h = fac * (m/rho)^(1/dim)
!           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
!           3 : other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )

    integer :: sle
    parameter(sle=0)

!     smoothing kernel function
!     skf = 1, cubic spline kernel by w4 - spline (monaghan 1985)
!         = 2, gauss kernel   (gingold and monaghan 1981)
!         = 3, quintic kernel (morris 1997)
    integer :: skf
    parameter(skf=1)

!     switches for different senarios

!     summation_density = .true. : use density summation model in the code,
!                        .false.: use continuiity equation
!     average_velocity = .true. : monaghan treatment on average velocity,
!                       .false.: no average treatment.
!     config_input = .true. : load initial configuration data,
!                   .false.: generate initial configuration.
!     virtual_part = .true. : use vritual particle,
!                   .false.: no use of vritual particle.
!     vp_input = .true. : load virtual particle information,
!               .false.: generate virtual particle information.
!     visc = .true. : consider viscosity,
!           .false.: no viscosity.
!     ex_force =.true. : consider external force,
!               .false.: no external force.
!     visc_artificial = .true. : consider artificial viscosity,
!                      .false.: no considering of artificial viscosity.
!     heat_artificial = .true. : consider artificial heating,
!                      .false.: no considering of artificial heating.
!     self_gravity = .true. : considering self_gravity,
!                    .false.: no considering of self_gravity
!     nor_density =  .true. : density normalization by using cspm,
!                    .false.: no normalization.
    logical :: summation_density, average_velocity, config_input, virtual_part, vp_input, visc, ex_force, &
        heat_artificial, visc_artificial, self_gravity, nor_density
    parameter(summation_density=.true.)
    parameter(average_velocity=.true.)
    parameter(config_input=.false.)
    parameter(virtual_part=.true.)
    parameter(vp_input=.false.)
    parameter(visc=.true.)
    parameter(ex_force=.true.)
    parameter(visc_artificial=.false.)
    parameter(heat_artificial=.false.)
    parameter(self_gravity=.false.)
    parameter(nor_density=.true.)

!     symmetry of the problem
!     nsym = 0 : no symmetry,
!          = 1 : axis symmetry,
!          = 2 : center symmetry.
    integer :: nsym
    parameter(nsym=0)

!     control parameters for output
!     int_stat = .true. : print statistics about sph particle interactions.
!                        including virtual particle information.
!     print_step: print timestep (on screen)
!     save_step : save timestep    (to disk file)
!     moni_particle: the particle number for information monitoring.
    logical :: int_stat
    parameter(int_stat=.true.)
    integer :: print_step, save_step, moni_particle
    parameter(print_step=100, save_step=500, moni_particle=1600)

    real(rk) :: pi
    parameter(pi=3.14159265358979323846)

!     simulation cases
!     shocktube = .true. : carry out shock tube simulation
!     shearcavity = .true. : carry out shear cavity simulation
    logical :: shocktube, shearcavity
    parameter(shocktube=.false.)
    parameter(shearcavity=.true.)
