!> including file for parameters and constants used
!> in the entire sph software packages.

module parameter

    use config_m, only: rk

    !@todo#2022-4: move dim/maxn/max_interaction to toml! 
    !> 求解问题的维度
    !> dim : dimension of the problem (1, 2 or 3)
    integer, parameter :: dim = 2

    !     parameters for the computational geometry,
    !     x_maxgeom : upper limit of allowed x-regime
    !     x_mingeom : lower limit of allowed x-regime
    !     y_maxgeom : upper limit of allowed y-regime
    !     y_mingeom : lower limit of allowed y-regime
    !     z_maxgeom : upper limit of allowed z-regime
    !     z_mingeom : lower limit of allowed z-regime
    real(rk) :: x_maxgeom, x_mingeom, y_maxgeom, y_mingeom, z_maxgeom, z_mingeom
    parameter(x_maxgeom=10._rk, x_mingeom=-10._rk, y_maxgeom=10._rk, y_mingeom=-10._rk, z_maxgeom=10._rk, z_mingeom=-10._rk)

    !> 粒子近似算法的指示变量
    !> sph algorithm for particle approximation (pa_sph)
    !> pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
    !>          2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
    integer, parameter :: pa_sph = 2

    !> 光滑长度估算方法的指示变量
    !> smoothing length evolution (sle) algorithm
    !> sle = 0 : keep unchanged,
    !>       1 : h = fac * (m/rho)^(1/dim)
    !>       2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
    !>       3 : other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
    integer, parameter :: sle = 0

    !     switches for different senarios

    !     summation_density = .true. : use density summation model in the code,
    !                        .false. : use continuiity equation
    !     average_velocity = .true. : monaghan treatment on average velocity,
    !                       .false. : no average treatment.
    !     config_input = .true. : load initial configuration data,
    !                   .false. : generate initial configuration.
    !     virtual_part = .true. : use vritual particle,
    !                   .false. : no use of vritual particle.
    !     vp_input = .true. : load virtual particle information,
    !               .false. : generate virtual particle information.
    !     visc = .true. : consider viscosity,
    !           .false. : no viscosity.
    !     ex_force =.true. : consider external force,
    !              .false. : no external force.
    !     visc_artificial = .true. : consider artificial viscosity,
    !                      .false. : no considering of artificial viscosity.
    !     heat_artificial = .true. : consider artificial heating,
    !                      .false. : no considering of artificial heating.
    
    ! 密度正则化的指示变量
    !     nor_density =  .true. : density normalization by using cspm,
    !                   .false. : no normalization.
    logical :: summation_density, average_velocity, config_input, virtual_part, vp_input, visc, ex_force, &
               heat_artificial, visc_artificial, nor_density
    parameter(summation_density=.true.)
    parameter(average_velocity=.true.)
    parameter(config_input=.false.)
    parameter(virtual_part=.true.)
    parameter(vp_input=.false.)
    parameter(visc=.true.)
    parameter(ex_force=.true.)
    parameter(visc_artificial=.false.)
    parameter(heat_artificial=.false.)
    parameter(nor_density=.true.)

    !     symmetry of the problem
    !     nsym = 0 : no symmetry,
    !          = 1 : axis symmetry,
    !          = 2 : center symmetry.
    integer :: nsym
    parameter(nsym=0)

    !> 控制 SPH 粒子相互作用状态的输出
    !> control parameters for output
    !> int_stat = .true. : print statistics about sph particle interactions.
    !>                    including virtual particle information.
    logical, parameter :: int_stat = .true.

    !> 所要监测的粒子的序号
    !> moni_particle: the particle number for information monitoring.
    integer, parameter :: moni_particle = 1600

    real(rk), parameter :: pi = acos(-1.0_rk)

    !> simulation cases

    !> 一维振荡管
    !> shocktube = .true. : carry out shock tube simulation
    logical, parameter :: shocktube = .false.

    !> 剪切腔
    !> shearcavity = .true. : carry out shear cavity simulation
    logical, parameter :: shearcavity = .true.

    !> 树型搜索边界尺度比例
    integer, parameter :: tree_scale_ratio = 10

end module parameter
