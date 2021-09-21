!>      x-- coordinates of particles                       [input/output]
!>      vx-- velocities of particles                       [input/output]
!>      mass-- mass of particles                                  [input]
!>      rho-- dnesities of particles                       [input/output]
!>      p-- pressure  of particles                         [input/output]
!>      u-- internal energy of particles                   [input/output]
!>      c-- sound velocity of particles                          [output]
!>      s-- entropy of particles, not used here                  [output]
!>      e-- total energy of particles                            [output]
!>      itype-- types of particles                               [input]
!>           =1   ideal gas
!>           =2   water
!>           =3   tnt
!>      hsml-- smoothing lengths of particles              [input/output]
!>      ntotal-- total particle number                            [input]
!>      maxtimestep-- maximum timesteps                           [input]
!>      dt-- timestep                                             [input]

subroutine time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: itype(maxn), ntotal, maxtimestep
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), dt
    integer  :: i, j, k, itimestep, d, current_ts, nstart
    real(rk) :: x_min(dim, maxn), v_min(dim, maxn), u_min(maxn), rho_min(maxn), dx(dim, maxn), dvx(dim, maxn), &
                du(maxn), drho(maxn), av(dim, maxn), ds(maxn), t(maxn), tdsdt(maxn)
    real(rk) :: time, temp_rho, temp_u

    !> 初始化时间：0.
    time = 0.0_rk

    do i = 1, ntotal
        do d = 1, dim
            av(d, i) = 0._rk
        end do
    end do

    do itimestep = nstart + 1, nstart + maxtimestep

        current_ts = current_ts + 1
        if (mod(itimestep, print_step) == 0) then
            write (*, *) '______________________________________________'
            write (*, *) '  current number of time step =', itimestep, '     current time=', time + dt
            write (*, *) '______________________________________________'
        end if

        !     if not first time step, then update thermal energy, density and
        !     velocity half a time step

        if (itimestep /= 1) then

            do i = 1, ntotal
                u_min(i) = u(i)
                temp_u = 0._rk
                if (dim == 1) temp_u = -nsym*p(i)*vx(1, i)/x(1, i)/rho(i)
                u(i) = u(i) + (dt/2.)*(du(i) + temp_u)
                if (u(i) < 0) u(i) = 0._rk

                if (.not. summation_density) then
                    rho_min(i) = rho(i)
                    temp_rho = 0._rk
                    if (dim == 1) temp_rho = -nsym*rho(i)*vx(1, i)/x(1, i)
                    rho(i) = rho(i) + (dt/2.)*(drho(i) + temp_rho)
                end if

                do d = 1, dim
                    v_min(d, i) = vx(d, i)
                       vx(d, i) = vx(d, i) + (dt/2.)*dvx(d, i)
                end do
            end do

        end if

        !---  definition of variables out of the function vector:

        call single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, tdsdt, dx, dvx, du, ds, drho, itype, av)

        if (itimestep == 1) then

            do i = 1, ntotal
                temp_u = 0._rk
                if (dim == 1) temp_u = -nsym*p(i)*vx(1, i)/x(1, i)/rho(i)
                u(i) = u(i) + (dt/2.)*(du(i) + temp_u)
                if (u(i) < 0) u(i) = 0._rk

                if (.not. summation_density) then
                    temp_rho = 0._rk
                    if (dim == 1) temp_rho = -nsym*rho(i)*vx(1, i)/x(1, i)
                    rho(i) = rho(i) + (dt/2.)*(drho(i) + temp_rho)
                end if

                do d = 1, dim
                    vx(d, i) = vx(d, i) + (dt/2.)*dvx(d, i) + av(d, i)
                     x(d, i) = x(d, i) + dt*vx(d, i)
                end do
            end do

        else

            do i = 1, ntotal
                temp_u = 0._rk
                if (dim == 1) temp_u = -nsym*p(i)*vx(1, i)/x(1, i)/rho(i)
                u(i) = u_min(i) + dt*(du(i) + temp_u)
                if (u(i) < 0) u(i) = 0._rk

                if (.not. summation_density) then
                    temp_rho = 0._rk
                    if (dim == 1) temp_rho = -nsym*rho(i)*vx(1, i)/x(1, i)
                    rho(i) = rho_min(i) + dt*(drho(i) + temp_rho)
                end if

                do d = 1, dim
                    vx(d, i) = v_min(d, i) + dt*dvx(d, i) + av(d, i)
                     x(d, i) = x(d, i) + dt*vx(d, i)
                end do
            end do

        end if

        time = time + dt

        if (mod(itimestep, save_step) == 0) then
            call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
        end if

        if (mod(itimestep, print_step) == 0) then
            write (*, *)
            write (*, 101) 'x', 'velocity', 'dvx'
            write (*, 100) x(1, moni_particle), vx(1, moni_particle), dvx(1, moni_particle)
        end if

    end do

    nstart = current_ts

101 format(1x, 3(2x, a12))
100 format(1x, 3(2x, e12.6))

end subroutine time_integration
