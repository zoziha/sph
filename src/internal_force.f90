!>   subroutine to calculate the internal forces on the right hand side
!>   of the navier-stokes equations, i.e. the pressure gradient and the
!>   gradient of the viscous stress tensor, used by the time integration.
!>   moreover the entropy production due to viscous dissipation, tds/dt,
!>   and the change of internal energy per mass, de/dt, are calculated.
!>
!>     itimestep: current timestep number                            [in]
!>     dt     :   time step                                          [in]
!>     ntotal : number of particles                                  [in]
!>     hsml   : smoothing length                                     [in]
!>     mass   : particle masses                                      [in]
!>     vx     : velocities of all particles                          [in]
!>     niac   : number of interaction pairs                          [in]
!>     rho    : density                                              [in]
!>     eta    : dynamic viscosity                                    [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     dwdx   : derivative of kernel with respect to x, y and z      [in]
!>     itype  : type of particle (material types)                    [in]
!>     u      : particle internal energy                             [in]
!>     x      : particle coordinates                                 [in]
!>     itype  : particle type                                        [in]
!>     t      : particle temperature                             [in/out]
!>     c      : particle sound speed                                [out]
!>     p      : particle pressure                                   [out]
!>     dvxdt  : acceleration with respect to x, y and z             [out]
!>     tdsdt  : production of viscous entropy                       [out]
!>     dedt   : change of specific internal energy                  [out]

subroutine int_force(itimestep, dt, ntotal, hsml, mass, vx, niac, rho, eta, pair_i, pair_j, &
                     dwdx, u, itype, x, t, c, p, dvxdt, tdsdt, dedt)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: itimestep, ntotal, niac, pair_i(max_interaction), pair_j(max_interaction), itype(maxn)
    real(rk) :: dt, hsml(maxn), mass(maxn), vx(dim, maxn), rho(maxn), eta(maxn), dwdx(dim, max_interaction), &
                u(maxn), x(dim, maxn), t(maxn), c(maxn), p(maxn), dvxdt(dim, maxn), tdsdt(maxn), dedt(maxn)
    integer  :: i, j, k, d
    real(rk) :: dvx(dim), txx(maxn), tyy(maxn), tzz(maxn), txy(maxn), txz(maxn), tyz(maxn), &
                vcc(maxn), hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij

    !     initialization of shear tensor, velocity divergence,
    !     viscous energy, internal energy, acceleration

    do i = 1, ntotal
        txx(i) = 0._rk
        tyy(i) = 0._rk
        tzz(i) = 0._rk
        txy(i) = 0._rk
        txz(i) = 0._rk
        tyz(i) = 0._rk
        vcc(i) = 0._rk
        tdsdt(i) = 0._rk
        dedt(i) = 0._rk
        do d = 1, dim
            dvxdt(d, i) = 0._rk
        end do
    end do

    !     calculate sph sum for shear tensor tab = va,b + vb,a - 2/3 delta_ab vc,c

    if (visc) then
        do k = 1, niac
            i = pair_i(k)
            j = pair_j(k)
            do d = 1, dim
                dvx(d) = vx(d, j) - vx(d, i)
            end do
            if (dim == 1) then
                hxx = 2._rk*dvx(1)*dwdx(1, k)
            else if (dim == 2) then
                hxx = 2._rk*dvx(1)*dwdx(1, k) - dvx(2)*dwdx(2, k)
                hxy = dvx(1)*dwdx(2, k) + dvx(2)*dwdx(1, k)
                hyy = 2._rk*dvx(2)*dwdx(2, k) - dvx(1)*dwdx(1, k)
            else if (dim == 3) then
                hxx = 2._rk*dvx(1)*dwdx(1, k) - dvx(2)*dwdx(2, k) - dvx(3)*dwdx(3, k)
                hxy = dvx(1)*dwdx(2, k) + dvx(2)*dwdx(1, k)
                hxz = dvx(1)*dwdx(3, k) + dvx(3)*dwdx(1, k)
                hyy = 2._rk*dvx(2)*dwdx(2, k) - dvx(1)*dwdx(1, k) - dvx(3)*dwdx(3, k)
                hyz = dvx(2)*dwdx(3, k) + dvx(3)*dwdx(2, k)
                hzz = 2._rk*dvx(3)*dwdx(3, k) - dvx(1)*dwdx(1, k) - dvx(2)*dwdx(2, k)
            end if
            hxx = 2._rk/3._rk*hxx
            hyy = 2._rk/3._rk*hyy
            hzz = 2._rk/3._rk*hzz
            if (dim == 1) then
                txx(i) = txx(i) + mass(j)*hxx/rho(j)
                txx(j) = txx(j) + mass(i)*hxx/rho(i)
            else if (dim == 2) then
                txx(i) = txx(i) + mass(j)*hxx/rho(j)
                txx(j) = txx(j) + mass(i)*hxx/rho(i)
                txy(i) = txy(i) + mass(j)*hxy/rho(j)
                txy(j) = txy(j) + mass(i)*hxy/rho(i)
                tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
                tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
            else if (dim == 3) then
                txx(i) = txx(i) + mass(j)*hxx/rho(j)
                txx(j) = txx(j) + mass(i)*hxx/rho(i)
                txy(i) = txy(i) + mass(j)*hxy/rho(j)
                txy(j) = txy(j) + mass(i)*hxy/rho(i)
                txz(i) = txz(i) + mass(j)*hxz/rho(j)
                txz(j) = txz(j) + mass(i)*hxz/rho(i)
                tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
                tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
                tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
                tyz(j) = tyz(j) + mass(i)*hyz/rho(i)
                tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
                tzz(j) = tzz(j) + mass(i)*hzz/rho(i)
            end if

            !     calculate sph sum for vc,c = dvx/dx + dvy/dy + dvz/dz:

            hvcc = 0.
            do d = 1, dim
                hvcc = hvcc + dvx(d)*dwdx(d, k)
            end do
            vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
            vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
        end do
    end if

    do i = 1, ntotal

        !     viscous entropy tds/dt = 1/2 eta/rho tab tab

        if (visc) then
            if (dim == 1) then
                tdsdt(i) = txx(i)*txx(i)
            else if (dim == 2) then
                tdsdt(i) = txx(i)*txx(i) + 2._rk*txy(i)*txy(i) + tyy(i)*tyy(i)
            else if (dim == 3) then
                tdsdt(i) = txx(i)*txx(i) + 2._rk*txy(i)*txy(i) + 2._rk*txz(i)*txz(i) + &
                           tyy(i)*tyy(i) + 2._rk*tyz(i)*tyz(i) + tzz(i)*tzz(i)
            end if
            tdsdt(i) = 0.5_rk*eta(i)/rho(i)*tdsdt(i)
        end if

        !     pressure from equation of state

        if (abs(itype(i)) == 1) then
            call p_gas(rho(i), u(i), p(i), c(i))
        else if (abs(itype(i)) == 2) then
            call p_art_water(rho(i), p(i), c(i))
        end if

    end do

    !      calculate sph sum for pressure force -p,a/rho
    !      and viscous force (eta tab),b/rho
    !      and the internal energy change de/dt due to -p/rho vc,c

    do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        he = 0._rk

        !     for sph algorithm 1

        rhoij = 1._rk/(rho(i)*rho(j))
        if (pa_sph == 1) then
            do d = 1, dim

                !     pressure part

                h = -(p(i) + p(j))*dwdx(d, k)
                he = he + (vx(d, j) - vx(d, i))*h

                !     viscous force

                if (visc) then

                    if (d == 1) then

                        !     x-coordinate of acceleration

                        h = h + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1, k)
                        if (dim >= 2) then
                            h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2, k)
                            if (dim == 3) then
                                h = h + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(3, k)
                            end if
                        end if
                    else if (d == 2) then

                        !     y-coordinate of acceleration

                        h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(1, k) + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2, k)
                        if (dim == 3) then
                            h = h + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(3, k)
                        end if
                    else if (d == 3) then

                        !     z-coordinate of acceleration

                        h = h + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(1, k) + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(2, k) &
                            + (eta(i)*tzz(i) + eta(j)*tzz(j))*dwdx(3, k)
                    end if
                end if
                h = h*rhoij
                dvxdt(d, i) = dvxdt(d, i) + mass(j)*h
                dvxdt(d, j) = dvxdt(d, j) - mass(i)*h
            end do
            he = he*rhoij
            dedt(i) = dedt(i) + mass(j)*he
            dedt(j) = dedt(j) + mass(i)*he

            !     for sph algorithm 2

        else if (pa_sph == 2) then
            do d = 1, dim
                h = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(d, k)
                he = he + (vx(d, j) - vx(d, i))*h

                !     viscous force

                if (visc) then
                    if (d == 1) then

                        !     x-coordinate of acceleration

                        h = h + (eta(i)*txx(i)/rho(i)**2 + eta(j)*txx(j)/rho(j)**2)*dwdx(1, k)
                        if (dim >= 2) then
                            h = h + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(2, k)
                            if (dim == 3) then
                                h = h + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(3, k)
                            end if
                        end if
                    else if (d == 2) then

                        !     y-coordinate of acceleration

                        h = h + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(1, k) + &
                            (eta(i)*tyy(i)/rho(i)**2 + eta(j)*tyy(j)/rho(j)**2)*dwdx(2, k)
                        if (dim == 3) then
                            h = h + (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(3, k)
                        end if
                    else if (d == 3) then

                        !     z-coordinate of acceleration

                        h = h + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(1, k) + &
                            (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(2, k) + &
                            (eta(i)*tzz(i)/rho(i)**2 + eta(j)*tzz(j)/rho(j)**2)*dwdx(3, k)
                    end if
                end if
                dvxdt(d, i) = dvxdt(d, i) + mass(j)*h
                dvxdt(d, j) = dvxdt(d, j) - mass(i)*h
            end do
            dedt(i) = dedt(i) + mass(j)*he
            dedt(j) = dedt(j) + mass(i)*he
        end if
    end do

    !     change of specific internal energy de/dt = t ds/dt - p/rho vc,c:

    do i = 1, ntotal
        dedt(i) = tdsdt(i) + 0.5_rk*dedt(i)
    end do

end subroutine int_force
