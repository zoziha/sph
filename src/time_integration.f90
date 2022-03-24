!> 时间步长前进，而且在整个系统中时间步长保持为常量。
!> 但是时间步长可以是时间和空间 (对应于每个粒子) 的变量。
!> 相关参考为 Hernquist 和 Katz (1989), Simpson (1995), Monaghan (1992) 等等。
subroutine time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt)

    use config_m, only: rk, stdout
    use parameter
    use output_m, only: output_all
    use progress_bar_m, only: pbflush, pbout
    use M_attr, only: fg_cyan, reset
    use console_color_m, only: attr, is_intel_windows
    implicit none

    !> 粒子的坐标
    !> coordinates of particles
    real(rk), intent(inout) :: x(dim, maxn)
    !> 粒子的速度
    !> velocities of particles
    real(rk), intent(inout) :: vx(dim, maxn)
    !> 粒子的质量
    !> mass of particles
    real(rk), intent(inout) :: mass(maxn)
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

    integer :: i, itimestep, d, nstart = 0  !! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量
    real(rk) :: x_min(dim, maxn), v_min(dim, maxn), u_min(maxn), rho_min(maxn), dx(dim, maxn), dvx(dim, maxn), &
                du(maxn), drho(maxn), ds(maxn), t(maxn), tdsdt(maxn)
    real(rk) :: av(dim, maxn)   !! 平均速度, average velocity
    real(rk) :: temp_rho, temp_u, &
                time = 0.0_rk   !! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量

    do itimestep = nstart + 1, nstart + maxtimestep    !! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量

        if (mod(itimestep, print_step) == 0) then
            call pbflush()      ! 进度条辅助程序
            write (stdout, '(a,i0)') attr('<INFO>')//'Current number of time step = ', itimestep
            write (stdout, "(a,g0.3,a)") attr('<INFO>')//'Current time = ', time, 's'
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

            !> 输出每个保存时间步的求解信息（拓展）
            call output_all(x, vx, mass, rho, p, u, c, itype, hsml, ntotal, itimestep/save_step)

        end if

        if (mod(itimestep, print_step) == 0) then
            write (*, *)
            write (*, 101) 'location', 'velocity', 'acc'
            write (*, 100) x(1, moni_particle), vx(1, moni_particle), dvx(1, moni_particle)
            !> 屏幕输出进度条
            if (.not.is_intel_windows()) write (stdout, '(a)', advance='no') fg_cyan
            call pbout(itimestep, nstart + maxtimestep, .true.)
            if (.not.is_intel_windows()) write (stdout, '(a)', advance='no') reset
            
        end if

    end do

    nstart = nstart + maxtimestep

101 format(3(a12,:,2x))
100 format(3(es12.5,:,2x))

end subroutine time_integration
