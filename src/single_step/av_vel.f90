!> 计算校正平均速度的子程序。详见 Monaghan(1992) (XSPH) 和第 4 章中的论述。
!>请问这在书中的多少页啊，找不到！希望能注释一下。
!> subroutine to calculate the average velocity to correct velocity
!> for preventing.penetration (Monaghan, 1992)
subroutine av_vel(ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)

    use sph_kind, only: rk
    use parameter
    implicit none

    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 粒子的质量
    !> particle masses
    real(rk) :: mass(maxn)
    !> 相互作用对的数目
    !> number of interaction pairs
    integer, intent(in) :: niac
    !> 相互作用对的第一部分的列表
    !> list of first partner of interaction pair
    integer, intent(in) :: pair_i(max_interaction)
    !> 相互作用对的第二部分的列表
    !> list of second partner of interaction pair
    integer, intent(in) :: pair_j(max_interaction)
    !> 给定相互作用对的光滑核函数
    !> kernel for all interaction pairs
    real(rk), intent(in) :: w(max_interaction)
    !> 粒子的速度
    !> particle velocities
    real(rk), intent(in) :: vx(dim, maxn)
    !> 粒子的密度
    !> particle densities
    real(rk), intent(in) :: rho(maxn)
    !> 粒子的平均速度
    !> average velocity of each particle
    real(rk), intent(out) :: av(dim, maxn)

    integer :: i, j, k, d
    real(rk) :: vcc, dvx(dim), epsilon
    
    ! epsilon -- 经验性的容差值，可能会导致不稳定
    ! 例如，一维震荡管问题，e (epsilon) <= 0.3
    !     epsilon --- a small constants chosen by experience, may lead to instability.
    !     for example, for the 1 dimensional shock tube problem, the e <= 0.3

    epsilon = 0.3_rk

    do i = 1, ntotal
        do d = 1, dim
            av(d, i) = 0.0_rk
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
