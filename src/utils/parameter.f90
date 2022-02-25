!> including file for parameters and constants used
!> in the entire sph software packages.

module parameter

    use sph_kinds, only: rk

    !> 求解问题的维度
    !> dim : dimension of the problem (1, 2 or 3)
    integer, parameter :: dim = 2

    !> 粒子的最大数目
    !> maxn    : maximum number of particles (gfortran: *<WARNING>* `-fmax-stack-var-size=`, r change the code to use an ALLOCATABLE array)
    integer, parameter :: maxn = 3000
    !> 粒子相互作用对的最大数目
    !> max_interation : maximum number of interaction pairs
    integer, parameter :: max_interaction = 20*maxn

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

    !> 最近相邻粒子搜索 (NNPS) 法的指示变量
    !> nearest neighbor particle searching (nnps) method
    !> nnps = 1 : simplest and direct searching
    !>        2 : sorting grid linked list
    !>        3 : tree algorithm
    integer, parameter :: nnps = 2

    !> 光滑长度估算方法的指示变量
    !> smoothing length evolution (sle) algorithm
    !> sle = 0 : keep unchanged,
    !>       1 : h = fac * (m/rho)^(1/dim)
    !>       2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
    !>       3 : other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
    integer, parameter :: sle = 0

    !> 光滑核函数的指示变量
    !> smoothing kernel function (skf)
    !> skf = 1, cubic spline kernel by w4 - spline (monaghan 1985)
    !>     = 2, gauss kernel   (gingold and monaghan 1981)
    !>     = 3, quintic kernel (morris 1997)
    integer, parameter :: skf = 1

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
    !     self_gravity = .true. : considering self_gravity,
    !                   .false. : no considering of self_gravity
    
    ! 密度正则化的指示变量
    !     nor_density =  .true. : density normalization by using cspm,
    !                   .false. : no normalization.
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

    !> 控制 SPH 粒子相互作用状态的输出
    !> control parameters for output
    !> int_stat = .true. : print statistics about sph particle interactions.
    !>                    including virtual particle information.
    logical, parameter :: int_stat = .true.

    !> 控制在屏幕上显示的粒子信息是第几个时间步
    !> print_step: print timestep (on screen)
    integer, parameter :: print_step = 100
    !> 控制保存到外部磁盘的粒子信息是第几个时间步
    !> save_step : save timestep  (to disk file)
    integer, parameter :: save_step = 500

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
