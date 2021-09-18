!>   subroutine to calculate the smoothing funciton for each particle and
!>   the interaction parameters used by the sph algorithm. interaction
!>   pairs are determined by directly comparing the particle distance
!>   with the corresponding smoothing length.
!>   see p.148 in chapter 4
!>
!>     itimestep : current time step                                 [in]
!>     ntotal    : number of particles                               [in]
!>     hsml      : smoothing length                                  [in]
!>     x         : coordinates of all particles                      [in]
!>     niac      : number of interaction pairs                      [out]
!>     pair_i    : list of first partner of interaction pair        [out]
!>     pair_j    : list of second partner of interaction pair       [out]
!>     w         : kernel for all interaction pairs                 [out]
!>     dwdx      : derivative of kernel with respect to x, y and z  [out]
!>     countiac  : number of neighboring particles                  [out]

subroutine direct_find(itimestep, ntotal, hsml, x, niac, pair_i, pair_j, w, dwdx, countiac)
   
    use sph_kind, only: rk
    implicit none
    include 'param.inc'

    integer  :: itimestep, ntotal, niac, pair_i(max_interaction), pair_j(max_interaction), countiac(maxn)
    real(rk) :: hsml(maxn), x(dim, maxn), w(max_interaction), dwdx(dim, max_interaction)
    integer  :: i, j, d, sumiac, maxiac, miniac, noiac, maxp, minp, scale_k
    real(rk) :: dxiac(dim), driac, r, mhsml, tdwdx(dim)
    !     smoothing kernel function
    !     skf = 1, cubic spline kernel by w4 - spline (monaghan 1985)
    !         = 2, gauss kernel   (gingold and monaghan 1981)
    !         = 3, quintic kernel (morris 1997)
    if      (skf == 1) then; scale_k = 2
    else if (skf == 2) then; scale_k = 3
    else if (skf == 3) then; scale_k = 3
    end if

    countiac(1:ntotal) = 0

    niac = 0

    do i = 1, ntotal - 1
        do j = i + 1, ntotal
            dxiac(1) = x(1, i) - x(1, j)
            driac = dxiac(1)*dxiac(1)
            do d = 2, dim
                dxiac(d) = x(d, i) - x(d, j)
                driac = driac + dxiac(d)*dxiac(d)
            end do
            mhsml = (hsml(i) + hsml(j))/2.
            if (sqrt(driac) < scale_k*mhsml) then
                if (niac < max_interaction) then

                    !     neighboring pair list, and totalinteraction number and
                    !     the interaction number for each particle

                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    r = sqrt(driac)
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1

                    !     kernel and derivations of kernel

                    call kernel(r, dxiac, mhsml, w(niac), tdwdx)
                    do d = 1, dim
                        dwdx(d, niac) = tdwdx(d)
                    end do
                else
                    error stop ' >>> error <<< : too many interactions'
                end if
            end if
        end do
    end do

    !     statistics for the interaction

    sumiac = 0
    maxiac = 0
    miniac = 1000
    noiac  = 0
    do i = 1, ntotal
        sumiac = sumiac + countiac(i)
        if (countiac(i) > maxiac) then
            maxiac = countiac(i)
            maxp = i
        end if
        if (countiac(i) < miniac) then
            miniac = countiac(i)
            minp = i
        end if
        if (countiac(i) == 0) noiac = noiac + 1
    end do

    if (mod(itimestep, print_step) == 0) then
        if (int_stat) then
            print *, ' >> statistics: interactions per particle:'
            print *, '**** particle:', maxp, ' maximal interactions:', maxiac
            print *, '**** particle:', minp, ' minimal interactions:', miniac
            print *, '**** average :', real(sumiac)/real(ntotal)
            print *, '**** total pairs : ', niac
            print *, '**** particles with no interactions:', noiac
        end if
    end if

end subroutine direct_find
