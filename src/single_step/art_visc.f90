!> 人工粘度
module art_visc_m

    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: art_visc

contains

    !> 计算人工粘度的子程序。详见 Monaghan (1992), Hernquist 和 Katz (1989) 或第 4 章中的论述。
    pure subroutine art_visc(ntotal, hsml, mass, x, vx, niac, rho, c, pair_i, pair_j, w, dwdx, dvxdt, dedt)
        integer, intent(in) :: ntotal       !! 在模拟中所使用的粒子总数
        real(rk), intent(in) :: hsml(:)     !! 光滑长度
        real(rk), intent(in) :: mass(:)     !! 粒子的质量
        real(rk), intent(in) :: x(:, :)     !! 粒子的坐标
        real(rk), intent(in) :: vx(:, :)    !! 粒子的速度
        integer, intent(in) :: niac         !! 相互作用对的数目
        real(rk), intent(in) :: rho(:)      !! 密度
        real(rk), intent(in) :: c(:)        !! 温度
        integer, intent(in) :: pair_i(:)    !! 相互作用对的第一个粒子
        integer, intent(in) :: pair_j(:)    !! 相互作用对的第二个粒子
        real(rk), intent(in) :: w(:)        !! 相互作用对的核函数
        real(rk), intent(in) :: dwdx(:, :)  !! 相互作用对的核函数的导数
        real(rk), intent(out) :: dvxdt(:, :)!! 相互作用对的加速度
        real(rk), intent(out) :: dedt(:)    !! 改变特定内部能量

        integer :: i, j, k, d
        real(rk) :: dx, dvx(dim), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml

        ! 人工粘度参数：
        ! 剪切粘度
        !     parameter for the artificial viscosity:
        !     shear viscosity
        parameter(alpha=1.0_rk)

        ! 散装粘度
        !     bulk viscosity
        parameter(beta=1.0_rk)

        ! 参数以避免奇点
        !     parameter to avoid singularities
        parameter(etq=0.1_rk)

        do i = 1, ntotal
            do d = 1, dim
                dvxdt(d, i) = 0.0_rk
            end do
            dedt(i) = 0.0_rk
        end do

        !     calculate sph sum for artificial viscosity

        do k = 1, niac
            i = pair_i(k)
            j = pair_j(k)
            mhsml = (hsml(i) + hsml(j))/2.0_rk
            vr = 0.0_rk
            rr = 0.0_rk
            do d = 1, dim
                dvx(d) = vx(d, i) - vx(d, j)
                dx = x(d, i) - x(d, j)
                vr = vr + dvx(d)*dx
                rr = rr + dx*dx
            end do

            !     artificial viscous force only if v_ij * r_ij < 0

            if (vr < 0.0_rk) then

                !     calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )

                muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)

                !     calculate piv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

                mc = 0.5_rk*(c(i) + c(j))
                mrho = 0.5_rk*(rho(i) + rho(j))
                piv = (beta*muv - alpha*mc)*muv/mrho

                !     calculate sph sum for artificial viscous force

                do d = 1, dim
                    h = -piv*dwdx(d, k)
                    dvxdt(d, i) = dvxdt(d, i) + mass(j)*h
                    dvxdt(d, j) = dvxdt(d, j) - mass(i)*h
                    dedt(i) = dedt(i) - mass(j)*dvx(d)*h
                    dedt(j) = dedt(j) - mass(i)*dvx(d)*h
                end do
            end if
        end do

        !     change of specific internal energy:

        do i = 1, ntotal
            dedt(i) = 0.5_rk*dedt(i)
        end do
    end subroutine art_visc

end module art_visc_m
