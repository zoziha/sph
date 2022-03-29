module art_heat_m
    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: art_heat
contains
    !> 计算人工热量的子程序。详见 Monaghan (1992), Fulk (1994) 或第 4 章中的论述。
    !> subroutine to calculate the artificial heat(fulk, 1994, p, a-17)
    !> see equ.(4.74)
    subroutine art_heat(ntotal, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w, dwdx, dedt)

        !> 在模拟中所使用的粒子总数
        !> number of particles in simulation
        integer, intent(in) :: ntotal
        !> 光滑长度
        !> smoothing length
        real(rk), intent(in) :: hsml(:)
        !> 粒子的质量
        !> particle masses
        real(rk), intent(in) :: mass(:)
        !> 粒子的坐标
        !> particle coordinates
        real(rk), intent(in) :: x(:, :)
        !> 粒子的速度
        !> particle velocities
        real(rk), intent(in) :: vx(:, :)
        !> 粒子的密度
        !> particle density
        real(rk), intent(in) :: rho(:)
        !> 粒子的特殊内部能量
        !> particle specific internal energy
        real(rk), intent(in) :: u(:)
        !> 粒子的声速
        !> particle sound velocity
        real(rk), intent(in) :: c(:)
        !> 粒子之间的互动对数
        !> number of interaction pairs
        integer, intent(in) :: niac
        !> 粒子之间的互动对的第一个粒子
        !> first partner of interaction pair
        integer, intent(in) :: pair_i(:)
        !> 粒子之间的互动对的第二个粒子
        !> second partner of interaction pair
        integer, intent(in) :: pair_j(:)
        !> 互动对的核函数
        !> kernel for all interaction pairs
        real(rk), intent(in) :: w(:)
        !> 互动对的核函数的导数
        !> derivative of kernel with respect to x, y and z
        real(rk), intent(in) :: dwdx(:, :)
        !> 生成的人工热量
        !> produced artificial heat, adding to energy eq.
        real(rk), intent(out) :: dedt(:)

        integer :: i, j, k, d
        real(rk) :: dx, dvx(dim), rr, h, mrho, mhsml, vcc(maxn), hvcc, mui, muj, muij, rdwdx, g1, g2

        !---  parameter for the artificial heat conduction:

        g1 = 0.1_rk
        g2 = 1.0_rk
        do i = 1, ntotal
            vcc(i) = 0.0_rk
            dedt(i) = 0.0_rk
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
            mhsml = (hsml(i) + hsml(j))/2.0_rk
            mrho = 0.5_rk*(rho(i) + rho(j))
            rr = 0.0_rk
            rdwdx = 0.0_rk
            do d = 1, dim
                dx = x(d, i) - x(d, j)
                rr = rr + dx*dx
                rdwdx = rdwdx + dx*dwdx(d, k)
            end do
            mui = g1*hsml(i)*c(i) + g2*hsml(i)**2*(abs(vcc(i)) - vcc(i))
            muj = g1*hsml(j)*c(j) + g2*hsml(j)**2*(abs(vcc(j)) - vcc(j))
            muij = 0.5_rk*(mui + muj)
            h = muij/(mrho*(rr + 0.01_rk*mhsml**2))*rdwdx
            dedt(i) = dedt(i) + mass(j)*h*(u(i) - u(j))
            dedt(j) = dedt(j) + mass(i)*h*(u(j) - u(i))

        end do

        do i = 1, ntotal
            dedt(i) = 2.0_rk*dedt(i)
        end do
    end subroutine art_heat

end module art_heat_m
