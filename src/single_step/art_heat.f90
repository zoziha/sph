!> 人工热量
module art_heat_m

    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: art_heat

contains

    !> 计算人工热量的子程序。详见 Monaghan (1992), Fulk (1994) 或第 4 章中的论述 (式 4.74)。
    pure subroutine art_heat(ntotal, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w, dwdx, dedt)
        integer, intent(in) :: ntotal       !! 粒子总数
        integer, intent(in) :: niac         !! 粒子对数量
        integer, intent(in) :: pair_i(:)    !! 粒子对的第一个粒子编号
        integer, intent(in) :: pair_j(:)    !! 粒子对的第二个粒子编号
        real(rk), intent(in) :: hsml(:)     !! 粒子光滑长度
        real(rk), intent(in) :: mass(:)     !! 粒子质量
        real(rk), intent(in) :: x(:, :)     !! 粒子位置
        real(rk), intent(in) :: vx(:, :)    !! 粒子速度
        real(rk), intent(in) :: rho(:)      !! 粒子密度
        real(rk), intent(in) :: u(:)        !! 粒子内能
        real(rk), intent(in) :: c(:)        !! 粒子速度的光束
        real(rk), intent(in) :: w(:)        !! 粒子核函数值
        real(rk), intent(in) :: dwdx(:, :)  !! 粒子核函数的导数
        real(rk), intent(out) :: dedt(:)    !! 粒子热量

        integer :: k
        real(rk) :: dvx(dim), rr, h, mrho, mhsml, vcc(ntotal), rdwdx
        real(rk), parameter :: g1 = 0.1_rk, g2 = 1.0_rk

        vcc = 0.0_rk
        dedt = 0.0_rk

        do k = 1, niac

            associate (i => pair_i(k), j => pair_j(k))
                dvx(:) = vx(:, j) - vx(:, i)
                associate (hvcc => sum(dvx(:)*dwdx(:, k)))
                    vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
                    vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
                end associate
            end associate

        end do

        do k = 1, niac

            associate (i => pair_i(k), j => pair_j(k))
                mhsml = (hsml(i) + hsml(j))/2
                mrho = (rho(i) + rho(j))/2
                associate (dx => x(:, i) - x(:, j))
                    rr = sum(dx*dx)
                    rdwdx = sum(dx*dwdx(:, k))
                end associate

                associate (mui => g1*hsml(i)*c(i) + g2*hsml(i)**2*(abs(vcc(i)) - vcc(i)), &
                           muj => g1*hsml(j)*c(j) + g2*hsml(j)**2*(abs(vcc(j)) - vcc(j)))
                    associate (h => (mui + muj)/2/(mrho*(rr + 0.01_rk*mhsml*mhsml))*rdwdx)
                        dedt(i) = dedt(i) + mass(j)*h*(u(i) - u(j))
                        dedt(j) = dedt(j) + mass(i)*h*(u(j) - u(i))
                    end associate
                end associate

            end associate

        end do

        dedt = 2*dedt
    end subroutine art_heat

end module art_heat_m
