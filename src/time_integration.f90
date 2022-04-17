!> 时间步迭代
module time_integration_m

    use art_heat_m, only: art_heat
    use art_visc_m, only: art_visc
    use av_vel_m, only: av_vel
    use config_m, only: rk, stdout, nnps, print_step, save_step, dofile, visc, &
                        visc_artificial, heat_artificial, virtual_part
    use density_m, only: sum_density, con_density
    use direct_find_m, only: direct_find
    use easy_string_m, only: to_string
    use external_force_m, only: ext_force
    use hsml_m, only: h_upgrade
    use info_m, only: operator(.c.)
    use input_m, only: virt_part
    use internal_force_m, only: int_force
    use output_m, only: output_all
    use link_list_m, only: link_list
    use lua_call_m, only: lua_virt_part
    use parameter
    use progress_bar_m, only: pbflush, pbout
    use stdlib_logger, only: stdlog => global_logger
    use tree_search_m, only: tree_search
    use viscosity_m, only: viscosity
    implicit none
    private

    public :: time_integration

contains

    !> 时间步长前进，而且在整个系统中时间步长保持为常量。
    !> 但是时间步长可以是时间和空间 (对应于每个粒子) 的变量。
    !> 相关参考为 Hernquist 和 Katz (1989), Simpson (1995), Monaghan (1992) 等等。
    subroutine time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt)
        use config_m, only: max_interaction
        real(rk), intent(inout) :: x(:, :)      !! 粒子的坐标
        real(rk), intent(inout) :: vx(:, :)     !! 粒子的速度
        real(rk), intent(inout) :: mass(:)      !! 粒子的质量
        real(rk), intent(inout) :: rho(:)       !! 粒子的密度
        real(rk), intent(inout) :: p(:)         !! 粒子的压力
        real(rk), intent(inout) :: u(:)         !! 粒子的内部能量
        real(rk), intent(out) :: c(:)           !! 粒子的声速
        real(rk), intent(out) :: s(:)           !! 粒子的熵
        real(rk), intent(out) :: e(:)           !! 粒子的总能量
        integer, intent(inout) :: itype(:)      !! 粒子的类型
        real(rk), intent(inout) :: hsml(:)      !! 粒子的平滑长度
        integer, intent(in) :: ntotal           !! 粒子的总数
        integer, intent(in) :: maxtimestep      !! 最大的时间步长
        real(rk), intent(in) :: dt              !! 时间步长

        integer :: i, itimestep, d, nvirt, nstart = 0  ! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量
        real(rk) :: x_min(dim, ntotal), v_min(dim, ntotal), u_min(ntotal), rho_min(ntotal), &
                    dx(dim, 2*ntotal), dvx(dim, 2*ntotal), &
                    du(2*ntotal), drho(2*ntotal), ds(2*ntotal), t(2*ntotal), tdsdt(2*ntotal)
        real(rk) :: av(dim, 2*ntotal)   ! 平均速度, average velocity
        real(rk) :: temp_rho, temp_u, &
                    time = 0.0_rk   ! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量

        do itimestep = nstart + 1, nstart + maxtimestep    !! 注意这里使用了Fortran的save属性，可以让程序在运行时保存这个变量

            if (mod(itimestep, print_step) == 0) then
                call pbflush()      ! 进度条辅助程序
                write (stdout, '(a,i0)') .c.'Current number of time step = ', itimestep
                write (stdout, "(a,g0.3,a)") .c.'Current time = ', time, 's'
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

            call single_step(max_interaction, itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, &
                             tdsdt, dx, dvx, du, ds, drho, itype, av, nvirt)

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

                !> 输出每个保存时间步的求解信息（拓展）, 输出实、虚粒子
                call output_all(x, vx, mass, rho, p, u, c, itype, hsml, ntotal + nvirt, itimestep/save_step)

            end if

            if (mod(itimestep, print_step) == 0) then
                write (*, *)
                write (*, 101) 'location', 'velocity', 'acc'
                do i = 1, dim
                    write (*, 100) x(i, moni_particle), vx(i, moni_particle), dvx(i, moni_particle)
                end do
                !> 屏幕输出进度条
                call pbout(itimestep, nstart + maxtimestep, isflushed=.true.)
            end if

        end do

        nstart = nstart + maxtimestep
        call stdlog%log_information('Total number of time steps, nstart = '//to_string(nstart))

101     format(3(a12, :, 2x))
100     format(3(es12.5, :, 2x))

    end subroutine time_integration

    !> 执行时间积分算法中的一个时间步的子程序
    subroutine single_step(max_interaction, itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, &
                           tdsdt, dx, dvx, du, ds, drho, itype, av, nvirt)
        integer, intent(in) :: max_interaction  !! 最大的互动数
        integer, intent(in) :: itimestep        !! 当前的时间步数
        real(rk), intent(in) :: dt              !! 当前的时间步长
        integer, intent(in) :: ntotal           !! 粒子总数
        real(rk), intent(inout) :: hsml(:)      !! 粒子的平滑长度
        real(rk), intent(inout) :: mass(:)      !! 粒子的质量
        real(rk), intent(inout) :: x(:, :)      !! 粒子的位置
        real(rk), intent(inout) :: vx(:, :)     !! 粒子的速度
        real(rk), intent(inout) :: u(:)         !! 粒子的内部能量
        real(rk), intent(inout) :: s(:)         !! 粒子的熵
        real(rk), intent(inout) :: rho(:)       !! 粒子的密度
        real(rk), intent(inout) :: p(:)         !! 粒子的压力
        real(rk), intent(inout) :: t(:)         !! 粒子的温度
        real(rk), intent(out) :: tdsdt(:)       !! 粒子的熵的生产量
        real(rk), intent(out) :: dx(:, :)       !! 粒子的位移的增量
        real(rk), intent(out) :: dvx(:, :)      !! 粒子的速度的增量
        real(rk), intent(out) :: du(:)          !! 粒子的内部能量的增量
        real(rk), intent(out) :: ds(:)          !! 粒子的熵的增量
        real(rk), intent(out) :: drho(:)        !! 粒子的密度的增量
        integer, intent(inout) :: itype(:)      !! 粒子的类型
        real(rk), intent(out) :: av(:, :)       !! 粒子的动量的增量
        integer, intent(out) :: nvirt           !! 粒子的虚粒子数量

        integer :: i
        logical, save :: loaded_virt = .false.
        !> 相互作用对的数目
        integer :: niac
        integer :: pair_i(max_interaction), pair_j(max_interaction), ns(2*ntotal) !@tofix: ntotal + nvirt, ifort, heap-array!
        real(rk) :: w(max_interaction), dwdx(dim, max_interaction), indvxdt(dim, 2*ntotal), &
                    exdvxdt(dim, 2*ntotal), ardvxdt(dim, 2*ntotal), avdudt(2*ntotal), ahdudt(2*ntotal), c(2*ntotal), eta(2*ntotal)

        do concurrent(i=1:ntotal)
            avdudt(i) = 0.0_rk
            ahdudt(i) = 0.0_rk
            indvxdt(:, i) = 0.0_rk
            ardvxdt(:, i) = 0.0_rk
        end do

        !> (边界) 虚粒子的位置设定
        !> positions of virtual (boundary) particles:

        nvirt = 0

        ! 初始化虚粒子的位置, 是否不需要重复设置虚粒子
        ! 从 Lua 脚本中获取设置虚粒子的方法或数据
        if (virtual_part) then
            if (dofile) then
                call lua_virt_part(x, vx, mass, rho, p, u, itype, hsml, ntotal, nvirt, keep=loaded_virt)
                if (.not. loaded_virt) loaded_virt = .true. !@tofix: more hack here needed
            else
                call virt_part(itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype, keep=loaded_virt)
                if (.not. loaded_virt) loaded_virt = .true.
            end if
        end if

        ! 交互作用参数，计算相邻粒子并优化平滑长度
        !---  interaction parameters, calculating neighboring particles
        !     and optimzing smoothing length

        if (nnps == 1) then
            call direct_find(itimestep, ntotal + nvirt, hsml, x, niac, pair_i, pair_j, w, dwdx, ns)
        elseif (nnps == 2) then
            call link_list(itimestep, ntotal + nvirt, hsml(1), x, niac, pair_i, pair_j, w, dwdx, ns)
        elseif (nnps == 3) then
            ! @todo: 测试1、3维算法
            call tree_search(itimestep, ntotal + nvirt, hsml, x, niac, pair_i, &
                             pair_j, w, dwdx, ns)
        end if

        ! 密度近似或改变rate，rate不知道如何翻译
        !---  density approximation or change rate

        if (summation_density) then
            call sum_density(ntotal + nvirt, hsml, mass, niac, pair_i, pair_j, w, itype, rho)
        else
            call con_density(ntotal + nvirt, mass, niac, pair_i, pair_j, dwdx, vx, itype, x, rho, drho)
        end if

        ! 动态粘性
        !---  dynamic viscosity:

        if (visc) call viscosity(ntotal + nvirt, itype, x, rho, eta)

        !---  internal forces:

        call int_force(itimestep, dt, ntotal + nvirt, hsml, mass, vx, niac, rho, eta, pair_i, pair_j, dwdx, &
                       u, itype, x, t, c, p, indvxdt, tdsdt, du)

        !---  artificial viscosity:

        if (visc_artificial) call art_visc(ntotal + nvirt, hsml, mass, x, vx, niac, rho, c, pair_i, pair_j, &
                                           w, dwdx, ardvxdt, avdudt)

        !---  external forces:

        if (ex_force) then
            call ext_force(ntotal + nvirt, mass, x, niac, pair_i, pair_j, itype, hsml, exdvxdt)
        else
            exdvxdt(:, 1:ntotal + nvirt) = 0.0_rk
        end if

        !     calculating the neighboring particles and undating hsml

        if (sle /= 0) call h_upgrade(dt, ntotal, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml)

        if (heat_artificial) call art_heat(ntotal + nvirt, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w, dwdx, ahdudt)

        ! 计算粒子的平均速度，避免粒子渗透 @todo: 重力作用下，仍会出现穿透
        !     calculating average velocity of each partile for avoiding penetration(渗透)

        if (average_velocity) then
            call av_vel(ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)
        else
            av(:, 1:ntotal) = 0.0_rk
        end if

        ! 转换粒子的速度，力和能量为f和dfdt
        !---  convert velocity, force, and energy to f and dfdt

        do concurrent(i=1:ntotal)
            dvx(:, i) = indvxdt(:, i) + exdvxdt(:, i) + ardvxdt(:, i)
            du(i) = du(i) + avdudt(i) + ahdudt(i)
        end do

        ! 监测粒子的第一个维度的速度改变量 (加速度)
        if (mod(itimestep, print_step) == 0) then
            write (stdout, 102) .c.'Information for particle, monitoring particle: ', moni_particle
            write (stdout, 101) 'internal a ', 'artifical a', 'external a =', 'total a '
            do i = 1, dim
                write (stdout, 100) indvxdt(i, moni_particle), ardvxdt(i, moni_particle), &
                    exdvxdt(i, moni_particle), dvx(i, moni_particle)
            end do
        end if

102     format(/a, i0)
101     format(4(a12, :, 2x))
100     format(4(es12.5, :, 2x))

    end subroutine single_step

end module time_integration_m
