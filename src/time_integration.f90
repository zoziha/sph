!> 时间步长前进，而且在整个系统中时间步长保持为常量。
!> 但是时间步长可以是时间和空间 (对应于每个粒子) 的变量。
!> 相关参考为 Hernquist 和 Katz (1989), Simpson (1995), Monaghan (1992) 等等。
subroutine time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt)

    use sph_kind, only: rk
    use parameter
    use output_m, only: output_all
    implicit none

    !> 粒子的坐标
    !> coordinates of particles
    real(rk), intent(inout) :: x(dim, maxn)
    !> 粒子的速度
    !> velocities of particles
    real(rk), intent(inout) :: vx(dim, maxn)
    !> 粒子的质量
    !> mass of particles
    real(rk), intent(in) :: mass(maxn)
    !> 粒子的密度
    !> dnesities of particles
    real(rk), intent(inout) :: rho(maxn)
    !> 粒子的压力
    !> pressure  of particles
    real(rk), intent(inout) :: p(maxn)
    !> 粒子的内部能量
    !> internal energy of particles
    real(rk), intent(inout) :: u(maxn)
    !> 粒子的声速
    !> sound velocity of particles
    real(rk), intent(out) :: c(maxn)
    !> 粒子的熵
    !> entropy of particles, not used here
    real(rk), intent(out) :: s(maxn)
    !> 粒子的总能量（@note: 暂未使用）
    !> total energy of particles
    real(rk), intent(out) :: e(maxn)
    !> 粒子的类型(1: ideal gas; 2: water; 3: TNT)
    !> types of particles
    integer, intent(inout) :: itype(maxn)
    !> 粒子的平滑长度
    !> smoothing lengths of particles
    real(rk), intent(inout) :: hsml(maxn)
    !> 粒子的总数
    !> total particle number
    integer, intent(in) :: ntotal
    !> 最大的时间步长
    !> maximum timesteps
    integer, intent(in) :: maxtimestep
    !> 时间步长
    !> timestep
    real(rk), intent(in) :: dt

    integer :: i, j, k, itimestep, d, nstart = 0  !! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量
    real(rk) :: x_min(dim, maxn), v_min(dim, maxn), u_min(maxn), rho_min(maxn), dx(dim, maxn), dvx(dim, maxn), &
                du(maxn), drho(maxn), ds(maxn), t(maxn), tdsdt(maxn)
    real(rk) :: av(dim, maxn)   !! 平均速度, average velocity
    real(rk) :: temp_rho, temp_u, &
                time = 0.0_rk   !! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量

    ! @todo: 重复初始化
    ! do i = 1, ntotal
    !     do d = 1, dim
    !         av(d, i) = 0.0_rk
    !     end do
    ! end do

    do itimestep = nstart + 1, nstart + maxtimestep    !! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量

        if (mod(itimestep, print_step) == 0) then
            write (*, *) '----------------------------------------------'
            write (*, '(1x,a,i0)') '  current number of time step = ', itimestep
            write (*, "(1x,a,g0.3)") '  current time = ', time
            write (*, *) '----------------------------------------------'
        end if

        ! 如果不是第一个时间步长，则更新热能、密度和速度半步长
        !     if not first time step, then update thermal energy, density and
        !     velocity half a time step

        if (itimestep /= 1) then

            do i = 1, ntotal
                u_min(i) = u(i)
                temp_u = 0.0_rk
                if (dim == 1) temp_u = -nsym*p(i)*vx(1, i)/x(1, i)/rho(i)
                u(i) = u(i) + (dt/2.0_rk)*(du(i) + temp_u)
                if (u(i) < 0) u(i) = 0.0_rk

                if (.not. summation_density) then
                    rho_min(i) = rho(i)
                    temp_rho = 0.0_rk
                    if (dim == 1) temp_rho = -nsym*rho(i)*vx(1, i)/x(1, i)
                    rho(i) = rho(i) + (dt/2.0_rk)*(drho(i) + temp_rho)
                end if

                do d = 1, dim
                    v_min(d, i) = vx(d, i)
                    vx(d, i) = vx(d, i) + (dt/2.0_rk)*dvx(d, i)
                end do
            end do

        end if

        !---  definition of variables out of the function vector:

        call single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, &
            tdsdt, dx, dvx, du, ds, drho, itype, av)

        if (itimestep == 1) then

            do i = 1, ntotal
                temp_u = 0.0_rk
                if (dim == 1) temp_u = -nsym*p(i)*vx(1, i)/x(1, i)/rho(i)
                u(i) = u(i) + (dt/2.0_rk)*(du(i) + temp_u)
                if (u(i) < 0) u(i) = 0.0_rk

                if (.not. summation_density) then
                    temp_rho = 0.0_rk
                    if (dim == 1) temp_rho = -nsym*rho(i)*vx(1, i)/x(1, i)
                    rho(i) = rho(i) + (dt/2.0_rk)*(drho(i) + temp_rho)
                end if

                ! 更新速度和位置
                do d = 1, dim
                    vx(d, i) = vx(d, i) + (dt/2.0_rk)*dvx(d, i) + av(d, i)
                    x(d, i) = x(d, i) + dt*vx(d, i)
                end do
            end do

        else

            do i = 1, ntotal
                temp_u = 0.0_rk
                if (dim == 1) temp_u = -nsym*p(i)*vx(1, i)/x(1, i)/rho(i)
                u(i) = u_min(i) + dt*(du(i) + temp_u)
                if (u(i) < 0) u(i) = 0.0_rk

                if (.not. summation_density) then
                    temp_rho = 0.0_rk
                    if (dim == 1) temp_rho = -nsym*rho(i)*vx(1, i)/x(1, i)
                    rho(i) = rho_min(i) + dt*(drho(i) + temp_rho)
                end if

                ! 更新速度和位置
                do d = 1, dim
                    vx(d, i) = v_min(d, i) + dt*dvx(d, i) + av(d, i)
                    x(d, i) = x(d, i) + dt*vx(d, i)
                end do
            end do

        end if

        time = time + dt

        if (mod(itimestep, save_step) == 0) then

            !> 覆盖输出最后的保存时间步的求解信息（局限）
            !> @todo: 待删除v2.0发布时
            !> @note: 重复覆盖输出最后一步的结果，意义不大，output_all足够了
            call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)

            !> 输出每个保存时间步的求解信息（拓展）
            call output_all(x, vx, mass, rho, p, u, c, itype, hsml, ntotal, itimestep/save_step)

        end if

        if (mod(itimestep, print_step) == 0) then
            write (*, *)
            write (*, 101) 'x', 'velocity', 'dvx'
            write (*, 100) x(1, moni_particle), vx(1, moni_particle), dvx(1, moni_particle)
        end if

    end do

    nstart = nstart + maxtimestep

101 format(1x, 3(2x, a12))
100 format(1x, 3(2x, es12.5))

end subroutine time_integration
