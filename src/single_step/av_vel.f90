!> 粒子速度校正
module av_vel_m

    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: av_vel

contains

    !> 计算校正平均速度的子程序。详见 Monaghan(1992) (XSPH) 和第 4 章中的论述。
    pure subroutine av_vel(ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)
        integer, intent(in) :: ntotal       !! 在模拟中所使用的粒子总数
        real(rk), intent(in) :: mass(:)     !! 粒子的质量
        integer, intent(in) :: niac         !! 相互作用对的数目
        integer, intent(in) :: pair_i(:)    !! 相互作用对的第一部分的列表
        integer, intent(in) :: pair_j(:)    !! 相互作用对的第二部分的列表
        real(rk), intent(in) :: w(:)        !! 给定相互作用对的光滑核函数
        real(rk), intent(in) :: vx(:, :)    !! 粒子的速度
        real(rk), intent(in) :: rho(:)      !! 粒子的密度
        real(rk), intent(out) :: av(:, :)   !! 粒子的平均速度

        integer :: i, j, k
        real(rk) :: dvx(dim), epsilon

        ! epsilon -- 经验性的容差值，可能会导致不稳定
        ! 例如，一维震荡管问题，e (epsilon) <= 0.3
        !     epsilon --- a small constants chosen by experience, may lead to instability.
        !     for example, for the 1 dimensional shock tube problem, the e <= 0.3
        parameter (epsilon = 0.3_rk)

        av(:, 1:ntotal) = 0.0_rk

        do k = 1, niac
            i = pair_i(k)
            j = pair_j(k)
            dvx = vx(:, i) - vx(:, j)
            av(:, i) = av(:, i) - 2*mass(j)*dvx(:)/(rho(i) + rho(j))*w(k)
            av(:, j) = av(:, j) + 2*mass(i)*dvx(:)/(rho(i) + rho(j))*w(k)
        end do

        av(:, 1:ntotal) = epsilon*av(:, 1:ntotal)
    end subroutine av_vel

end module av_vel_m
