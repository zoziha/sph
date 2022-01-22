!> 通过简单的粒子间距与光滑长度相匹配进行最近相邻粒子搜索的子程序。详见第 4 章 148 页。
!>   subroutine to calculate the smoothing funciton for each particle and
!>   the interaction parameters used by the sph algorithm. interaction
!>   pairs are determined by directly comparing the particle distance
!>   with the corresponding smoothing length.
!>   see p.148 in chapter 4
subroutine direct_find(itimestep, ntotal, hsml, x, niac, pair_i, pair_j, w, dwdx, countiac)

    use sph_kind, only: rk
    use parameter
    implicit none

    !> 当前时间步
    !> current time step
    integer, intent(in) :: itimestep
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 粒子的光滑长度
    !> smoothing length
    real(rk), intent(in) :: hsml(maxn)
    !> 粒子的坐标
    !> coordinates of all particles
    real(rk), intent(in) :: x(dim, maxn)
    !> 相互作用对的数目
    !> number of interaction pairs
    integer, intent(out) :: niac
    !> 相互作用对的第一个粒子
    !> first partner of interaction pair
    integer, intent(out) :: pair_i(max_interaction)
    !> 相互作用对的第二个粒子
    !> second partner of interaction pair
    integer, intent(out) :: pair_j(max_interaction)
    !> 给定相互作用对的光滑核函数
    !> kernel for all interaction pairs
    real(rk), intent(out) :: w(max_interaction)
    !> 核函数对 x, y, z 的导数
    !> derivative of kernel with respect to x, y and z
    real(rk), intent(out) :: dwdx(dim, max_interaction)
    !> 相互作用对的数目
    !> number of neighboring particles
    integer, intent(out) :: countiac(ntotal)

    integer :: i, j, d, sumiac, maxiac, miniac, noiac, maxp, minp, scale_k
    real(rk) :: dxiac(dim), driac, r, mhsml, tdwdx(dim)

    !> smoothing kernel function
    !> skf = 1, cubic spline kernel by w4 - spline (monaghan 1985)
    !>     = 2, gauss kernel   (gingold and monaghan 1981)
    !>     = 3, quintic kernel (morris 1997)
    if (skf == 1) then; scale_k = 2
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
            mhsml = (hsml(i) + hsml(j))/2.0_rk
            if (sqrt(driac) < scale_k*mhsml) then
                if (niac < max_interaction) then

                    !> neighboring pair list, and totalinteraction number and
                    !> the interaction number for each particle

                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    r = sqrt(driac)
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1

                    !> kernel and derivations of kernel

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

    !> statistics for the interaction

    sumiac = 0
    maxiac = 0
    miniac = 1000
    noiac = 0
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
            print 100, '**** particle: ', maxp, '   maximal interactions: ', maxiac
            print 100, '**** particle: ', minp, '   minimal interactions: ', miniac
            print 101, '**** average : ', real(sumiac)/real(ntotal)
            print 100, '**** total pairs : ', niac
            print 100, '**** particles with no interactions: ', noiac
        end if
    end if
    
101 format(1x,*(a,g0.2))
100 format(1x,*(a,i0))

end subroutine direct_find
