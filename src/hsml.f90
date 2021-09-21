!> Subroutine to evolve smoothing length.
!>
!>     dt     : time step                                            [in]
!>     ntotal : number of particles                                  [in]
!>     mass   : particle masses                                      [in]
!>     vx     : velocities of all particles                          [in]
!>     rho    : density                                              [in]
!>     niac   : number of interaction pairs                          [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     dwdx   : derivative of kernel with respect to x, y and z      [in]
!>     hsml   : smoothing length                                 [in/out]

subroutine h_upgrade(dt, ntotal, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: ntotal, niac, pair_i(max_interaction), pair_j(max_interaction)
    real(rk) :: mass(maxn), vx(dim, maxn), rho(maxn), dwdx(dim, max_interaction), hsml(maxn)
    integer  :: i, j, k, d
    real(rk) :: dt, fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn)

    if (sle == 0) then

        !---  keep smoothing length unchanged.

        return

    else if (sle == 2) then

        !---  dh/dt = (-1/dim)*(h/rho)*(drho/dt).

        do i = 1, ntotal
            vcc(i) = 0._rk
        end do

        do k = 1, niac
            i = pair_i(k)
            j = pair_j(k)
            do d = 1, dim
                dvx(d) = vx(d, j) - vx(d, i)
            end do
            hvcc = dvx(1)*dwdx(1, k)
            do d = 2, dim
                hvcc = hvcc + dvx(d)*dwdx(d, k)
            end do
            vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
            vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
        end do

        do i = 1, ntotal
            dhsml(i) = (hsml(i)/dim)*vcc(i)
             hsml(i) = hsml(i) + dt*dhsml(i)
            if (hsml(i) <= 0) hsml(i) = hsml(i) - dt*dhsml(i)
        end do

    else if (sle == 1) then

        fac = 2.0_rk
        do i = 1, ntotal
            hsml(i) = fac*(mass(i)/rho(i))**(1._rk/dim)
        end do

    end if

end subroutine h_upgrade
