!> 计算人工粘度的子程序。详见 Monaghan (1992), Hernquist 和 Katz (1989) 或第 4 章中的论述。
!>     subroutine to calculate the artificial viscosity (monaghan, 1992)
subroutine art_visc(ntotal, hsml, mass, x, vx, niac, rho, c, pair_i, pair_j, w, dwdx, dvxdt, dedt)

    use sph_kind, only: rk
    use parameter
    implicit none

    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 光滑长度
    !> smoothing length
    real(rk), intent(in) :: hsml(maxn)
    !> 粒子的质量
    !> particle masses
    real(rk), intent(in) :: mass(maxn)
    !> 粒子的坐标
    !> particle coordinates
    real(rk), intent(in) :: x(dim, maxn)
    !> 粒子的速度
    !> particle velocities
    real(rk), intent(in) :: vx(dim, maxn)
    !> 相互作用对的数目
    !> number of interaction pairs
    integer, intent(in) :: niac
    !> 密度
    !> density
    real(rk), intent(in) :: rho(maxn)

    !> 温度 ?
    !> temperature ?
    !> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!教材中说c是声速!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk), intent(in) :: c(maxn)
    !> !-----------------------------------------------------------------!

    !> 相互作用对的第一个粒子
    !> first partner of interaction pair
    integer, intent(in) :: pair_i(max_interaction)
    !> 相互作用对的第二个粒子
    !> second partner of interaction pair
    integer, intent(in) :: pair_j(max_interaction)
    !> 相互作用对的核函数
    !> kernel for all interaction pairs
    real(rk), intent(in) :: w(max_interaction)
    !> 相互作用对的核函数的导数
    !> derivative of kernel with respect to x, y and z
    real(rk), intent(in) :: dwdx(dim, max_interaction)
    !> 相互作用对的加速度
    !> acceleration with respect to x, y and z
    real(rk), intent(out) :: dvxdt(dim, maxn)
    !> 改变特定内部能量(我认为应该是能量改变率)
    !> change of specific internal energy
    real(rk), intent(out) :: dedt(maxn)

    integer :: i, j, k, d
    real(rk) :: dx, dvx(dim), alpha, beta, etq, piv, muv, vr, rr, h, mc, mrho, mhsml

    !     parameter for the artificial viscosity:
    !     shear viscosity
    parameter(alpha=1._rk)

    !     bulk viscosity
    parameter(beta=1._rk)

    !     parameter to avoid singularities
    parameter(etq=0.1_rk)

    do i = 1, ntotal
        do d = 1, dim
            dvxdt(d, i) = 0._rk
        end do
        dedt(i) = 0._rk
    end do

    !     calculate sph sum for artificial viscosity

    do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        mhsml = (hsml(i) + hsml(j))/2._rk !公式(4.70)
        vr = 0._rk
        rr = 0._rk
        do d = 1, dim
            dvx(d) = vx(d, i) - vx(d, j) !公式(4.71)
            dx = x(d, i) - x(d, j) !公式(4.71)
            vr = vr + dvx(d)*dx
            rr = rr + dx*dx
        end do

        !     artificial viscous force only if v_ij * r_ij < 0
        !只有当速度差分或距离差分小于零时才出现
        if (vr < 0._rk) then

            !     calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )

            muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)

            !     calculate piv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

            mc = 0.5_rk*(c(i) + c(j)) !公式(4.68)
            mrho = 0.5_rk*(rho(i) + rho(j)) !公式(4.69)
            piv = (beta*muv - alpha*mc)*muv/mrho !公式(4.66)，计算人工粘度

            !     calculate sph sum for artificial viscous force
            !公式(4.66)，计算人工粘度

            do d = 1, dim
                h = -piv*dwdx(d, k)
                dvxdt(d, i) = dvxdt(d, i) + mass(j)*h
                dvxdt(d, j) = dvxdt(d, j) - mass(i)*h !作用力与反作用力！
                dedt(i) = dedt(i) - mass(j)*dvx(d)*h  !这里的负号是因为冲击导致热量产生，导致能量耗散，故为负号
                dedt(j) = dedt(j) - mass(i)*dvx(d)*h  !同理 
            end do
        end if
    end do

    !     change of specific internal energy:

    do i = 1, ntotal
        dedt(i) = 0.5_rk*dedt(i)
    end do

end subroutine art_visc
