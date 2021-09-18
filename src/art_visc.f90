!>     subroutine to calculate the artificial viscosity (monaghan, 1992)
!>     see equ.(4.66) equ.(4.62)
!>
!>     ntotal : number of particles (including virtual particles)    [in]
!>     hsml   : smoothing length                                     [in]
!>     mass   : particle masses                                      [in]
!>     x      : coordinates of all particles                         [in]
!>     vx     : velocities of all particles                          [in]
!>     niac   : number of interaction pairs                          [in]
!>     rho    : density                                              [in]
!>     c      : temperature                                          [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     w      : kernel for all interaction pairs                     [in]
!>     dwdx   : derivative of kernel with respect to x, y and z      [in]
!>     dvxdt  : acceleration with respect to x, y and z             [out]
!>     dedt   : change of specific internal energy                  [out]

subroutine art_visc(ntotal, hsml, mass, x, vx, niac, rho, c, pair_i, pair_j, w, dwdx, dvxdt, dedt)

    use sph_kind, only: rk
    implicit none
    include 'param.inc'

    integer  :: ntotal, niac, pair_i(max_interaction), pair_j(max_interaction)
    real(rk) :: hsml(maxn), mass(maxn), x(dim, maxn), vx(dim, maxn), rho(maxn), c(maxn), w(max_interaction), &
                dwdx(dim, max_interaction), dvxdt(dim, maxn), dedt(maxn)
    integer  :: i, j, k, d
    real(rk) :: dx, dvx(dim), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml

    !     parameter for the artificial viscosity:
    !     shear viscosity
    parameter(alpha=1._rk)

    !     bulk viscosity
    parameter(beta=1._rk)

    !     parameter to avoid singularities
    parameter(etq =0.1_rk)

    do i = 1, ntotal
        do d = 1, dim
            dvxdt(d, i) = 0._rk
        end do
        dedt(i) = 0._rk
    end do

    !     calculate sph sum for artificial viscosity

    do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        mhsml = (hsml(i) + hsml(j))/2.
        vr = 0._rk
        rr = 0._rk
        do d = 1, dim
            dvx(d) = vx(d, i) - vx(d, j)
            dx = x(d, i) - x(d, j)
            vr = vr + dvx(d)*dx
            rr = rr + dx*dx
        end do

        !     artificial viscous force only if v_ij * r_ij < 0

        if (vr < 0._rk) then

            !     calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )

            muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)

            !     calculate piv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

            mc = 0.5_rk*(c(i) + c(j))
            mrho = 0.5_rk*(rho(i) + rho(j))
            piv = (beta*muv - alpha*mc)*muv/mrho

            !     calculate sph sum for artificial viscous force

            do d = 1, dim
                h = -piv*dwdx(d, k)
                dvxdt(d, i) = dvxdt(d, i) + mass(j)*h
                dvxdt(d, j) = dvxdt(d, j) - mass(i)*h
                dedt(i) = dedt(i) - mass(j)*dvx(d)*h
                dedt(j) = dedt(j) - mass(i)*dvx(d)*h
            end do
        end if
    end do

    !     change of specific internal energy:

    do i = 1, ntotal
        dedt(i) = 0.5_rk*dedt(i)
    end do

end subroutine art_visc
