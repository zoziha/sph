!> Subroutine to determine the right hand side of a differential
!>  equation in a single step for performing time integration.
!>
!> In this routine and its subroutines the sph algorithms are performed.
!>     itimestep: current timestep number                            [in]
!>     dt       : timestep                                           [in]
!>     ntotal   :  number of particles                               [in]
!>     hsml     :  smoothing length                                  [in]
!>     mass     :  particle masses                                   [in]
!>     x        :  particle position                                 [in]
!>     vx       :  particle velocity                                 [in]
!>     u        :  particle internal energy                          [in]
!>     s        :  particle entropy (not used here)                  [in]
!>     rho      :  density                                       [in/out]
!>     p        :  pressure                                         [out]
!>     t        :  temperature                                   [in/out]
!>     tdsdt    :  production of viscous entropy t*ds/dt            [out]
!>     dx       :  dx = vx = dx/dt                                  [out]
!>     dvx      :  dvx = dvx/dt, force per unit mass                [out]
!>     du       :  du  = du/dt                                      [out]
!>     ds       :  ds  = ds/dt                                      [out]
!>     drho     :  drho =  drho/dt                                  [out]
!>     itype    :  type of particle                                 [in]
!>     av       :  monaghan average velocity                        [out]

subroutine single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, tdsdt, dx, dvx, du, ds, drho, itype, av)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: itimestep, ntotal, itype(maxn)
    real(rk) :: dt, hsml(maxn), mass(maxn), x(dim, maxn), vx(dim, maxn), u(maxn), &
                s(maxn), rho(maxn), p(maxn), t(maxn), tdsdt(maxn), dx(dim, maxn), dvx(dim, maxn), &
                du(maxn), ds(maxn), drho(maxn), av(dim, maxn)
    integer  :: i, d, nvirt, niac, pair_i(max_interaction), pair_j(max_interaction), ns(maxn)
    real(rk) :: w(max_interaction), dwdx(dim, max_interaction), indvxdt(dim, maxn), &
                exdvxdt(dim, maxn), ardvxdt(dim, maxn), avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn)

    do i = 1, ntotal
        avdudt(i) = 0._rk
        ahdudt(i) = 0._rk
        do d = 1, dim
            indvxdt(d, i) = 0._rk
            ardvxdt(d, i) = 0._rk
            exdvxdt(d, i) = 0._rk
        end do
    end do

    !---  positions of virtual (boundary) particles:
    !> （边界）虚粒子的位置设定

    nvirt = 0
    if (virtual_part) then; call virt_part(itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype)
    end if

    !---  interaction parameters, calculating neighboring particles
    !     and optimzing smoothing length

    if      (nnps == 1) then; call direct_find(itimestep, ntotal + nvirt, hsml, x, niac, pair_i, pair_j, w, dwdx, ns)
    else if (nnps == 2) then; call link_list(itimestep, ntotal + nvirt, hsml(1), x, niac, pair_i, pair_j, w, dwdx, ns)
    else if (nnps == 3) then
        !        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
        !     &       pair_j,w,dwdx,ns)
    end if

    !---  density approximation or change rate

    if (summation_density) then; call sum_density(ntotal + nvirt, hsml, mass, niac, pair_i, pair_j, w, itype, rho)
    else;                        call con_density(ntotal + nvirt, mass, niac, pair_i, pair_j, dwdx, vx, itype, x, rho, drho)
    end if

    !---  dynamic viscosity:

    if (visc) call viscosity(ntotal + nvirt, itype, x, rho, eta)

    !---  internal forces:

    call int_force(itimestep, dt, ntotal + nvirt, hsml, mass, vx, niac, rho, eta, pair_i, pair_j, dwdx, &
                   u, itype, x, t, c, p, indvxdt, tdsdt, du)

    !---  artificial viscosity:

    if (visc_artificial) call art_visc(ntotal + nvirt, hsml, mass, x, vx, niac, rho, c, pair_i, pair_j, &
                                       w, dwdx, ardvxdt, avdudt)

    !---  external forces:

    if (ex_force) call ext_force(ntotal + nvirt, mass, x, niac, pair_i, pair_j, itype, hsml, exdvxdt)

    !     calculating the neighboring particles and undating hsml

    if (sle /= 0) call h_upgrade(dt, ntotal, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml)

    if (heat_artificial) call art_heat(ntotal + nvirt, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w, dwdx, ahdudt)

    !     calculating average velocity of each partile for avoiding penetration

    if (average_velocity) call av_vel(ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)

    !---  convert velocity, force, and energy to f and dfdt

    do i = 1, ntotal
        do d = 1, dim
            dvx(d, i) = indvxdt(d, i) + exdvxdt(d, i) + ardvxdt(d, i)
        end do
        du(i) = du(i) + avdudt(i) + ahdudt(i)
    end do

    if (mod(itimestep, print_step) == 0) then
        write (*, *)
        write (*, *) '**** information for particle ****', moni_particle
        write (*, 101) 'internal a ', 'artifical a=', 'external a ', 'total a '
        write (*, 100) indvxdt(1, moni_particle), ardvxdt(1, moni_particle), exdvxdt(1, moni_particle), dvx(1, moni_particle)
    end if

101 format(1x, 4(2x, a12))
100 format(1x, 4(2x, e12.6))

end subroutine single_step
