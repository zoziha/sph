!>     subroutine to calculate the artificial heat(fulk, 1994, p, a-17)
!>     see equ.(4.74)
!>
!>     ntotal : number of particles                                  [in]
!>     hsml   : smoothing length                                     [in]
!>     mass   : particle masses                                      [in]
!>     x      : coordinates of all particles                         [in]
!>     vx     : velocities of all particles                          [in]
!>     rho    : density                                              [in]
!>     u      : specific internal energy                             [in]
!>     c      : sound veolcity                                       [in]
!>     niac   : number of interaction pairs                          [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     w      : kernel for all interaction pairs                     [in]
!>     dwdx   : derivative of kernel with respect to x, y and z      [in]
!>     dedt   : produced artificial heat, adding to energy eq.      [out]

subroutine art_heat(ntotal, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w, dwdx, dedt)

    use sph_kind, only: rk
    implicit none
    include 'param.inc'

    integer  :: ntotal, niac, pair_i(max_interaction), pair_j(max_interaction)
    real(rk) :: hsml(maxn), mass(maxn), x(dim, maxn), vx(dim, maxn), rho(maxn), &
                u(maxn), c(maxn), w(max_interaction), dwdx(dim, max_interaction), dedt(maxn)
    integer  :: i, j, k, d
    real(rk) :: dx, dvx(dim), vr, rr, h, mc, mrho, mhsml, vcc(maxn), hvcc, mui, muj, muij, rdwdx, g1, g2

    !---  parameter for the artificial heat conduction:

    g1 = 0.1_rk
    g2 = 1.0_rk
    do i = 1, ntotal
        vcc(i)  = 0._rk
        dedt(i) = 0._rk
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

    do k = 1, niac

        i = pair_i(k)
        j = pair_j(k)
        mhsml = (hsml(i) + hsml(j))/2._rk
        mrho  = 0.5_rk*(rho(i) + rho(j))
        rr    = 0._rk
        rdwdx = 0._rk
        do d = 1, dim
            dx = x(d, i) - x(d, j)
            rr = rr + dx*dx
            rdwdx = rdwdx + dx*dwdx(d, k)
        end do
        mui  = g1*hsml(i)*c(i) + g2*hsml(i)**2*(abs(vcc(i)) - vcc(i))
        muj  = g1*hsml(j)*c(j) + g2*hsml(j)**2*(abs(vcc(j)) - vcc(j))
        muij = 0.5*(mui + muj)
        h = muij/(mrho*(rr + 0.01*mhsml**2))*rdwdx
        dedt(i) = dedt(i) + mass(j)*h*(u(i) - u(j))
        dedt(j) = dedt(j) + mass(i)*h*(u(j) - u(i))

    end do

    do i = 1, ntotal
        dedt(i) = 2.0_rk*dedt(i)
    end do

end subroutine art_heat
