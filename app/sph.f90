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

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: ntotal, itype(maxn), maxtimestep, d, m, i, yesorno
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), &
                u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), dt
    real(rk) :: s1, s2

    call time_print
    call cpu_time(s1)

    if (shocktube)   dt = 0.005_rk
    if (shearcavity) dt = 5.e-5_rk
    call input(x, vx, mass, rho, p, u, itype, hsml, ntotal)
1   write (*, *) '  ***************************************************'
    write (*, *) '          please input the maximal time steps '
    write (*, *) '  ***************************************************'
    read  (*, *) maxtimestep
    call time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt)

    !> 输出最后一个时间布的求解信息
    call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
    
    write (*, *) '  ***************************************************'
    write (*, *) 'are you going to run more time steps ? (0=no, 1=yes)'
    write (*, *) '  ***************************************************'
    read  (*, *) yesorno
    if (yesorno /= 0) goto 1
    call time_print
    call cpu_time(s2)
    write (*, *) '        elapsed cpu time (seconds) = ', s2 - s1
    write (*, *) 'all finish!'

end program sph
