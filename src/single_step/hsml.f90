!> 更新光滑长度
!> Subroutine to evolve smoothing length.
subroutine h_upgrade(dt, ntotal, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml)

    use sph_kind, only: rk
    use parameter
    implicit none

    !> 时间步长
    !> Time step
    real(rk), intent(in) :: dt
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 粒子的质量
    !> particle masses
    real(rk), intent(in) :: mass(maxn)
    !> 粒子的速度
    !> particle velocities
    real(rk), intent(in) :: vx(dim, maxn)
    !> 密度
    !> density
    real(rk), intent(in) :: rho(maxn)
    !> 相互作用对的数目
    !> number of interaction pairs
    integer, intent(in) :: niac
    !> 相互作用对的第一个粒子
    !> first partner of interaction pair
    integer, intent(in) :: pair_i(max_interaction)
    !> 相互作用对的第二个粒子
    !> second partner of interaction pair
    integer, intent(in) :: pair_j(max_interaction)
    !> 对应于粒子的每个方向的核函数导数
    !> derivative of kernel with respect to x, y and z
    real(rk), intent(in) :: dwdx(dim, max_interaction)
    !> 光滑长度
    !> smoothing length
    real(rk), intent(inout) :: hsml(maxn)

    integer :: i, j, k, d
    real(rk) :: fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn)

    !> 光滑长度估算方法的指示变量
    !> smoothing length evolution (sle) algorithm
    !> sle = 0 : keep unchanged,
    !>       1 : h = fac * (m/rho)^(1/dim)
    !>       2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
    !>       3 : other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
    if (sle == 0) then

        !---  keep smoothing length unchanged.

        return

    else if (sle == 2) then !公式（4.81）

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

    else if (sle == 1) then !公式（4.80）

        fac = 2.0_rk
        do i = 1, ntotal
            hsml(i) = fac*(mass(i)/rho(i))**(1._rk/dim)
        end do

    end if

end subroutine h_upgrade
