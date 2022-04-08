module internal_force_m

    use config_m, only: rk, visc
    use eos_m, only: p_gas, p_art_water
    use parameter, only: dim, pa_sph
    implicit none
    private
    
    public :: int_force

contains

    !> 计算SPH相互作用力的子程序。
    !> 内力的计算是通过剪切应力的嵌套 SPH 近似法进行的。
    !> 在此子程序中实现了两种类型的 SPH 粒子近似法。
    !> 详见第 4 章，其他相关参考有 Riffert 等人 (1995), Flebbe 等人 (1994)。
    !> Subroutine to calculate the internal forces on the right hand side
    !>  of the navier-stokes equations, i.e. the pressure gradient and the
    !>  gradient of the viscous stress tensor, used by the time integration.
    !>  moreover the entropy production due to viscous dissipation, tds/dt,
    !>  and the change of internal energy per mass, de/dt, are calculated.
    subroutine int_force(itimestep, dt, ntotal, hsml, mass, vx, niac, rho, eta, pair_i, pair_j, &
                         dwdx, u, itype, x, t, c, p, dvxdt, tdsdt, dedt)
        !> 当前时间步
        !> Current time step
        integer, intent(in) :: itimestep
        !> 时间步长
        !> Time step
        real(rk), intent(in) :: dt
        !> 粒子数量
        !> Number of particles
        integer, intent(in) :: ntotal
        !> 光滑长度
        !> Smoothing length
        real(rk), intent(in) :: hsml(:)
        !> 粒子质量
        !> Particle masses
        real(rk), intent(in) :: mass(:)
        !> 粒子速度
        !> Particle velocities
        real(rk), intent(in) :: vx(:, :)
        !> 粒子互动对数
        !> Number of interaction pairs
        integer, intent(in) :: niac
        !> 粒子密度
        !> Particle density
        real(rk), intent(in) :: rho(:)
        !> 动态粘性
        !> Dynamic viscosity
        real(rk), intent(in) :: eta(:)
        !> 互动对第一个粒子
        !> First partner of interaction pair
        integer, intent(in) :: pair_i(:)
        !> 互动对第二个粒子
        !> Second partner of interaction pair
        integer, intent(in) :: pair_j(:)
        !> 核函数对于x, y, z的导数
        !> Derivative of kernel with respect to x, y and z
        real(rk), intent(in) :: dwdx(:, :)
        !> 粒子内部能量
        !> Particle internal energy
        real(rk), intent(in) :: u(:)
        !> 粒子坐标
        !> Particle coordinates
        real(rk), intent(in) :: x(:, :)
        !> 粒子类型
        !> Particle type
        integer, intent(in) :: itype(:)
        !> 粒子温度
        !> Particle temperature
        real(rk), intent(inout) :: t(:)
        !> 粒子声速
        !> Particle sound speed
        real(rk), intent(out) :: c(:)
        !> 粒子压力
        !> Particle pressure
        real(rk), intent(out) :: p(:)
        !> 粒子加速度
        !> Acceleration with respect to x, y and z
        real(rk), intent(out) :: dvxdt(:, :)
        !> 粒子消耗的粘性熵
        !> Production of viscous entropy
        real(rk), intent(out) :: tdsdt(:)
        !> 粒子温度变化
        !> Change of specific internal energy
        real(rk), intent(out) :: dedt(:)

        integer :: i, j, k, d
        real(rk) :: dvx(dim), txx(ntotal), tyy(ntotal), tzz(ntotal), txy(ntotal), txz(ntotal), tyz(ntotal), &
                    vcc(ntotal), hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij

        !     initialization of shear tensor, velocity divergence,
        !     viscous energy, internal energy, acceleration

        do i = 1, ntotal
            txx(i) = 0._rk
            tyy(i) = 0._rk
            tzz(i) = 0._rk
            txy(i) = 0._rk
            txz(i) = 0._rk
            tyz(i) = 0._rk
            vcc(i) = 0._rk
            tdsdt(i) = 0._rk
            dedt(i) = 0._rk
            do d = 1, dim
                dvxdt(d, i) = 0._rk
            end do
        end do

        !     calculate sph sum for shear tensor tab = va,b + vb,a - 2/3 delta_ab vc,c

        if (visc) then
            do k = 1, niac
                i = pair_i(k)
                j = pair_j(k)
                do d = 1, dim
                    dvx(d) = vx(d, j) - vx(d, i)
                end do
                if (dim == 1) then
                    hxx = 2._rk*dvx(1)*dwdx(1, k)
                else if (dim == 2) then
                    hxx = 2._rk*dvx(1)*dwdx(1, k) - dvx(2)*dwdx(2, k)
                    hxy = dvx(1)*dwdx(2, k) + dvx(2)*dwdx(1, k)
                    hyy = 2._rk*dvx(2)*dwdx(2, k) - dvx(1)*dwdx(1, k)
                else if (dim == 3) then
                    hxx = 2._rk*dvx(1)*dwdx(1, k) - dvx(2)*dwdx(2, k) - dvx(3)*dwdx(3, k)
                    hxy = dvx(1)*dwdx(2, k) + dvx(2)*dwdx(1, k)
                    hxz = dvx(1)*dwdx(3, k) + dvx(3)*dwdx(1, k)
                    hyy = 2._rk*dvx(2)*dwdx(2, k) - dvx(1)*dwdx(1, k) - dvx(3)*dwdx(3, k)
                    hyz = dvx(2)*dwdx(3, k) + dvx(3)*dwdx(2, k)
                    hzz = 2._rk*dvx(3)*dwdx(3, k) - dvx(1)*dwdx(1, k) - dvx(2)*dwdx(2, k)
                end if
                hxx = 2._rk/3._rk*hxx
                hyy = 2._rk/3._rk*hyy
                hzz = 2._rk/3._rk*hzz
                if (dim == 1) then
                    txx(i) = txx(i) + mass(j)*hxx/rho(j)
                    txx(j) = txx(j) + mass(i)*hxx/rho(i)
                else if (dim == 2) then
                    txx(i) = txx(i) + mass(j)*hxx/rho(j)
                    txx(j) = txx(j) + mass(i)*hxx/rho(i)
                    txy(i) = txy(i) + mass(j)*hxy/rho(j)
                    txy(j) = txy(j) + mass(i)*hxy/rho(i)
                    tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
                    tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
                else if (dim == 3) then
                    txx(i) = txx(i) + mass(j)*hxx/rho(j)
                    txx(j) = txx(j) + mass(i)*hxx/rho(i)
                    txy(i) = txy(i) + mass(j)*hxy/rho(j)
                    txy(j) = txy(j) + mass(i)*hxy/rho(i)
                    txz(i) = txz(i) + mass(j)*hxz/rho(j)
                    txz(j) = txz(j) + mass(i)*hxz/rho(i)
                    tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
                    tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
                    tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
                    tyz(j) = tyz(j) + mass(i)*hyz/rho(i)
                    tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
                    tzz(j) = tzz(j) + mass(i)*hzz/rho(i)
                end if

                !     calculate sph sum for vc,c = dvx/dx + dvy/dy + dvz/dz:

                hvcc = 0._rk
                do d = 1, dim
                    hvcc = hvcc + dvx(d)*dwdx(d, k)
                end do
                vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
                vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
            end do
        end if

        do i = 1, ntotal

            !     viscous entropy tds/dt = 1/2 eta/rho tab tab

            if (visc) then
                if (dim == 1) then
                    tdsdt(i) = txx(i)*txx(i)
                else if (dim == 2) then
                    tdsdt(i) = txx(i)*txx(i) + 2._rk*txy(i)*txy(i) + tyy(i)*tyy(i)
                else if (dim == 3) then
                    tdsdt(i) = txx(i)*txx(i) + 2._rk*txy(i)*txy(i) + 2._rk*txz(i)*txz(i) + &
                               tyy(i)*tyy(i) + 2._rk*tyz(i)*tyz(i) + tzz(i)*tzz(i)
                end if
                tdsdt(i) = 0.5_rk*eta(i)/rho(i)*tdsdt(i)
            end if

            !     pressure from equation of state

            if (abs(itype(i)) == 1) then
                call p_gas(rho(i), u(i), p(i), c(i))
            else if (abs(itype(i)) == 2) then
                call p_art_water(rho(i), p(i), c(i))
            end if

        end do

        !      calculate sph sum for pressure force -p,a/rho
        !      and viscous force (eta tab),b/rho
        !      and the internal energy change de/dt due to -p/rho vc,c

        do k = 1, niac
            i = pair_i(k)
            j = pair_j(k)
            he = 0._rk

            !     for sph algorithm 1

            rhoij = 1._rk/(rho(i)*rho(j))
            if (pa_sph == 1) then
                do d = 1, dim

                    !     pressure part

                    h = -(p(i) + p(j))*dwdx(d, k)
                    he = he + (vx(d, j) - vx(d, i))*h

                    !     viscous force

                    if (visc) then

                        if (d == 1) then

                            !     x-coordinate of acceleration

                            h = h + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1, k)
                            if (dim >= 2) then
                                h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2, k)
                                if (dim == 3) then
                                    h = h + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(3, k)
                                end if
                            end if
                        else if (d == 2) then

                            !     y-coordinate of acceleration

                            h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(1, k) + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2, k)
                            if (dim == 3) then
                                h = h + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(3, k)
                            end if
                        else if (d == 3) then

                            !     z-coordinate of acceleration

                            h = h + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(1, k) + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(2, k) &
                                + (eta(i)*tzz(i) + eta(j)*tzz(j))*dwdx(3, k)
                        end if
                    end if
                    h = h*rhoij
                    dvxdt(d, i) = dvxdt(d, i) + mass(j)*h
                    dvxdt(d, j) = dvxdt(d, j) - mass(i)*h
                end do
                he = he*rhoij
                dedt(i) = dedt(i) + mass(j)*he
                dedt(j) = dedt(j) + mass(i)*he

                !     for sph algorithm 2

            else if (pa_sph == 2) then
                do d = 1, dim
                    h = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(d, k)
                    he = he + (vx(d, j) - vx(d, i))*h

                    !     viscous force

                    if (visc) then
                        if (d == 1) then

                            !     x-coordinate of acceleration

                            h = h + (eta(i)*txx(i)/rho(i)**2 + eta(j)*txx(j)/rho(j)**2)*dwdx(1, k)
                            if (dim >= 2) then
                                h = h + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(2, k)
                                if (dim == 3) then
                                    h = h + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(3, k)
                                end if
                            end if
                        else if (d == 2) then

                            !     y-coordinate of acceleration

                            h = h + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(1, k) + &
                                (eta(i)*tyy(i)/rho(i)**2 + eta(j)*tyy(j)/rho(j)**2)*dwdx(2, k)
                            if (dim == 3) then
                                h = h + (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(3, k)
                            end if
                        else if (d == 3) then

                            !     z-coordinate of acceleration

                            h = h + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(1, k) + &
                                (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(2, k) + &
                                (eta(i)*tzz(i)/rho(i)**2 + eta(j)*tzz(j)/rho(j)**2)*dwdx(3, k)
                        end if
                    end if
                    dvxdt(d, i) = dvxdt(d, i) + mass(j)*h
                    dvxdt(d, j) = dvxdt(d, j) - mass(i)*h
                end do
                dedt(i) = dedt(i) + mass(j)*he
                dedt(j) = dedt(j) + mass(i)*he
            end if
        end do

        !     change of specific internal energy de/dt = t ds/dt - p/rho vc,c:

        do i = 1, ntotal
            dedt(i) = tdsdt(i) + 0.5_rk*dedt(i)
        end do
    end subroutine int_force

end module internal_force_m
