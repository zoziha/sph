!>     subroutine to calculate the average velocity to correct velocity
!>     for preventing.penetration (monaghan, 1992)
!>
!>     ntotal : number of particles                                  [in]
!>     mass   : particle masses                                      [in]
!>     niac   : number of interaction pairs                          [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     w      : kernel for all interaction pairs                     [in]
!>     vx     : velocity of each particle                            [in]
!>     rho    : density of each particle                             [in]
!>     av     : average velocityof each particle                    [out]

subroutine av_vel(ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: ntotal, niac, pair_i(max_interaction), pair_j(max_interaction)
    real(rk) :: mass(maxn), w(max_interaction), vx(dim, maxn), rho(maxn), av(dim, maxn)
    integer  :: i, j, k, d
    real(rk) :: vcc, dvx(dim), epsilon

    !     epsilon --- a small constants chosen by experience, may lead to instability.
    !     for example, for the 1 dimensional shock tube problem, the e <= 0.3

    epsilon = 0.3

    do i = 1, ntotal
        do d = 1, dim
            av(d, i) = 0.
        end do
    end do

    do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        do d = 1, dim
            dvx(d) = vx(d, i) - vx(d, j)
            av(d, i) = av(d, i) - 2*mass(j)*dvx(d)/(rho(i) + rho(j))*w(k)
            av(d, j) = av(d, j) + 2*mass(i)*dvx(d)/(rho(i) + rho(j))*w(k)
        end do
    end do

    do i = 1, ntotal
        do d = 1, dim
            av(d, i) = epsilon*av(d, i)
        end do
    end do

end subroutine av_vel
