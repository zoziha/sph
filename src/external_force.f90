!>     subroutine to calculate the external forces, e.g. gravitational forces.
!>     the forces from the interactions with boundary virtual particles
!>     are also calculated here as external forces.
!>
!>     here as the external force.
!>     ntotal  : number of particles                                 [in]
!>     mass    : particle masses                                     [in]
!>     x       : coordinates of all particles                        [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     itype   : type of particles                                   [in]
!>     hsml   : smoothing length                                     [in]
!>     dvxdt   : acceleration with respect to x, y and z            [out]

subroutine ext_force(ntotal, mass, x, niac, pair_i, pair_j, itype, hsml, dvxdt)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: ntotal, itype(maxn), niac, pair_i(max_interaction), pair_j(max_interaction)
    real(rk) :: mass(maxn), x(dim, maxn), hsml(maxn), dvxdt(dim, maxn)
    integer  :: i, j, k, d
    real(rk) :: dx(dim), rr, f, rr0, dd, p1, p2

    do i = 1, ntotal
        do d = 1, dim
            dvxdt(d, i) = 0._rk
        end do
    end do

    !     consider self-gravity or not ?

    if (self_gravity) then
        do i = 1, ntotal
            dvxdt(dim, i) = -9.8_rk
        end do
    end if

    !     boundary particle force and penalty anti-penetration force.
    rr0 = 1.25e-5
    dd = 1.e-2
    p1 = 12
    p2 = 4

    do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        if (itype(i) > 0 .and. itype(j) < 0) then
            rr = 0.
            do d = 1, dim
                dx(d) = x(d, i) - x(d, j)
                rr = rr + dx(d)*dx(d)
            end do
            rr = sqrt(rr)
            if (rr < rr0) then
                f = ((rr0/rr)**p1 - (rr0/rr)**p2)/rr**2
                do d = 1, dim
                    dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
                end do
            end if
        end if
    end do

end subroutine ext_force
