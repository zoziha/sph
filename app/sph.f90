!> This is a three dimensional sph code. the followings are the
!>  basic parameters needed in this code or calculated by this code.
!>
!>     mass-- mass of particles                                      [in]
!>     ntotal-- total particle number ues                            [in]
!>     dt--- time step used in the time integration                  [in]
!>     itype-- types of particles                                    [in]
!>     x-- coordinates of particles                              [in/out]
!>     vx-- velocities of particles                              [in/out]
!>     rho-- densities of particles                              [in/out]
!>     p-- pressure  of particles                                [in/out]
!>     u-- internal energy of particles                          [in/out]
!>     hsml-- smoothing lengths of particles                     [in/out]
!>     c-- sound velocity of particles                              [out]
!>     s-- entropy of particles                                     [out]
!>     e-- total energy of particles                                [out]

program sph

    use sph_kinds, only: rk
    use parameter
    use utils, only: tic, toc
    use output_m, only: set_parameter_log
    implicit none

    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer :: ntotal
    !> 粒子的类型(1: ideal gas; 2: water)
    !> types of particles
    integer :: itype(maxn)
    integer :: maxtimestep
    integer :: d, m, i, yesorno
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), &
                u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), dt

    call tic()
    call time_print()
    call set_parameter_log()

    if (shocktube) dt = 0.005_rk
    if (shearcavity) dt = 5.0e-5_rk
    call input(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    ! 主循环
    do
        write (*, "(a)", advance="no") 'Please input the maximal time steps: '
        read (*, *) maxtimestep

        call time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt)

        !> 输出最后一个时间步的求解信息
        call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)

        write (*, "(a)", advance="no") 'Are you going to run more time steps ? (0=no, 1=yes): '
        read (*, *) yesorno
        if (yesorno == 0) exit

    end do
    call time_print()
    call toc()
    write (*, "(a)") 'All finish!'

end program sph
