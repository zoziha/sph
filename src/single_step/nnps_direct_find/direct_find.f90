module direct_find_m

    use config_m, only: rk, skf
    use parameter
    use utils, only: get_distance
    use output_m, only: set_statistics_print
    use kernel_m, only: kernel
    implicit none
    private

    public :: direct_find

contains

    !> 通过简单的粒子间距与光滑长度相匹配进行最近相邻粒子搜索的子程序。详见第 4 章 148 页。
    !>   subroutine to calculate the smoothing funciton for each particle and
    !>   the interaction parameters used by the sph algorithm. interaction
    !>   pairs are determined by directly comparing the particle distance
    !>   with the corresponding smoothing length.
    !>   see p.148 in chapter 4
    subroutine direct_find(itimestep, ntotal, hsml, x, niac, pair_i, pair_j, w, dwdx, countiac)

        !> 当前时间步
        !> current time step
        integer, intent(in) :: itimestep
        !> 在模拟中所使用的粒子总数
        !> number of particles in simulation
        integer, intent(in) :: ntotal
        !> 粒子的光滑长度
        !> smoothing length
        real(rk), intent(in) :: hsml(:)
        !> 粒子的坐标
        !> coordinates of all particles
        real(rk), intent(in) :: x(:, :)
        !> 相互作用对的数目
        !> number of interaction pairs
        integer, intent(out) :: niac
        !> 相互作用对的第一个粒子
        !> first partner of interaction pair
        integer, intent(out) :: pair_i(:)
        !> 相互作用对的第二个粒子
        !> second partner of interaction pair
        integer, intent(out) :: pair_j(:)
        !> 给定相互作用对的光滑核函数
        !> kernel for all interaction pairs
        real(rk), intent(out) :: w(:)
        !> 核函数对 x, y, z 的导数
        !> derivative of kernel with respect to x, y and z
        real(rk), intent(out) :: dwdx(:, :)
        !> 相互作用对的数目
        !> number of neighboring particles
        integer, intent(out) :: countiac(:)

        integer :: i, j, scale_k
        real(rk) :: dxiac(dim), r, mhsml

        !> 光滑核函数
        !> smoothing kernel function
        !> skf = 1, cubic spline kernel by w4 - spline (monaghan 1985)
        !>     = 2, gauss kernel   (gingold and monaghan 1981)
        !>     = 3, quintic kernel (morris 1997)
        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        countiac(1:ntotal) = 0
        niac = 0

        do i = 1, ntotal - 1
            do j = i + 1, ntotal

                ! 计算两个粒子之间的距离的平方
                ! calculate distance between two particles
                call get_distance(x(1:dim, i), x(1:dim, j), dxiac, r)

                ! 对称光滑长度: 光滑长度的算数平均值 Page 127.
                mhsml = (hsml(i) + hsml(j))/2
                if (r < scale_k*mhsml) then
                    if (niac < max_interaction) then

                        !> 相邻对列表，以及每个粒子的总交互次数和相互作用数
                        !> neighboring pair list, and totalinteraction number and
                        !> the interaction number for each particle

                        niac = niac + 1
                        pair_i(niac) = i
                        pair_j(niac) = j
                        countiac(i) = countiac(i) + 1
                        countiac(j) = countiac(j) + 1

                        !> 核函数及其对 x, y, z 的导数
                        !> kernel and derivations of kernel
                        call kernel(r, dxiac, mhsml, w(niac), dwdx(:, niac))

                    else
                        error stop ' >>> error <<< : too many interactions'
                    end if
                end if
            end do
        end do

        !> 相互作用的统计信息
        !     statistics for the interaction
        call set_statistics_print(itimestep, ntotal, niac, countiac)
    end subroutine direct_find
    
end module direct_find_m
