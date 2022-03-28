!> 执行时间积分算法中的一个时间步的子程序
!> Subroutine to determine the right hand side of a differential
!>  equation in a single step for performing time integration.
!>
!> In this routine and its subroutines the sph algorithms are performed.
subroutine single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, &
                       tdsdt, dx, dvx, du, ds, drho, itype, av)

    use config_m, only: rk, stdout, nnps
    use parameter
    use tree_search_m, only: tree_search
    use info_m, only: operator(.c.)
    implicit none

    !> 当前时间步
    !> Current timestep
    integer, intent(in) :: itimestep
    !> 时间步长
    !> Time step
    real(rk), intent(in) :: dt
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 粒子的平滑长度
    !> smoothing length
    real(rk), intent(inout) :: hsml(maxn)
    !> 粒子的质量
    !> particle masses
    real(rk), intent(inout) :: mass(maxn)
    !> 粒子的位置
    !> particle positions
    real(rk), intent(inout) :: x(dim, maxn)
    !> 粒子的速度
    !> particle velocities
    real(rk), intent(inout) :: vx(dim, maxn)
    !> 粒子的内部能量
    !> particle internal energy
    real(rk), intent(inout) :: u(maxn)
    !> 粒子的熵
    !> particle entropy (not used here)
    real(rk), intent(in) :: s(maxn)
    !> 粒子的密度
    !> particle density
    real(rk), intent(inout) :: rho(maxn)
    !> 粒子的压力
    !> particle pressure
    real(rk), intent(inout) :: p(maxn)
    !> 粒子的温度
    !> particle temperature
    real(rk), intent(inout) :: t(maxn)
    !> 粒子的熵的生产量
    !> production of viscous entropy t*ds/dt
    real(rk), intent(out) :: tdsdt(maxn)
    !> 粒子的位移的增量
    !> particle displacement, dx = vx = dx/dt
    real(rk), intent(out) :: dx(dim, maxn)
    !> 粒子的速度的增量
    !> particle velocity displacement, dvx = dvx/dt
    real(rk), intent(out) :: dvx(dim, maxn)
    !> 粒子的内部能量的增量
    !> particle internal energy displacement, du = du/dt
    real(rk), intent(out) :: du(maxn)
    !> 粒子的熵的增量
    !> particle entropy displacement, ds  = ds/dt
    real(rk), intent(out) :: ds(maxn)
    !> 粒子的密度的增量
    !> particle density displacement, drho =  drho/dt
    real(rk), intent(out) :: drho(maxn)
    !> 粒子的类型 (1: ideal gas; 2: water)
    !> particle type
    integer, intent(inout) :: itype(maxn)
    !> 平均速度
    !> monaghan average velocity
    real(rk), intent(out) :: av(dim, maxn)

    integer :: i, d, nvirt
    !> 相互作用对的数目
    integer :: niac
    integer :: pair_i(max_interaction), pair_j(max_interaction), ns(maxn)
    real(rk) :: w(max_interaction), dwdx(dim, max_interaction), indvxdt(dim, maxn), &
                exdvxdt(dim, maxn), ardvxdt(dim, maxn), avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn)

    do i = 1, ntotal
        avdudt(i) = 0.0_rk
        ahdudt(i) = 0.0_rk
        do d = 1, dim
            indvxdt(d, i) = 0.0_rk
            ardvxdt(d, i) = 0.0_rk
            exdvxdt(d, i) = 0.0_rk
        end do
    end do

    !> (边界) 虚粒子的位置设定
    !> positions of virtual (boundary) particles:

    nvirt = 0
    if (virtual_part) then
        call virt_part(itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype)
    end if

    ! 交互作用参数，计算相邻粒子并优化平滑长度
    !---  interaction parameters, calculating neighboring particles
    !     and optimzing smoothing length

    if (nnps == 1) then
        call direct_find(itimestep, ntotal + nvirt, hsml, x, niac, pair_i, pair_j, w, dwdx, ns)
    elseif (nnps == 2) then
        call link_list(itimestep, ntotal + nvirt, hsml(1), x, niac, pair_i, pair_j, w, dwdx, ns)
    elseif (nnps == 3) then
        ! @todo: 树型搜索算法（zoziha/quad-tree: https://github.com/zoziha/quad-tree）
        ! @tocheck
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

    if (ex_force) call ext_force(ntotal + nvirt, mass, x, niac, pair_i, pair_j, itype, hsml, exdvxdt)

    !     calculating the neighboring particles and undating hsml

    if (sle /= 0) call h_upgrade(dt, ntotal, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml)

    if (heat_artificial) call art_heat(ntotal + nvirt, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w, dwdx, ahdudt)

    ! 计算粒子的平均速度，避免粒子渗透
    !     calculating average velocity of each partile for avoiding penetration(渗透)

    if (average_velocity) call av_vel(ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)

    ! 转换粒子的速度，力和能量为f和dfdt
    !---  convert velocity, force, and energy to f and dfdt

    do i = 1, ntotal
        do d = 1, dim
            dvx(d, i) = indvxdt(d, i) + exdvxdt(d, i) + ardvxdt(d, i)
        end do
        du(i) = du(i) + avdudt(i) + ahdudt(i)
    end do

    ! 监测粒子的第一个维度的速度改变量 (加速度)
    if (mod(itimestep, print_step) == 0) then
        write (stdout, 102) .c.'Information for particle, monitoring particle: ', moni_particle
        write (stdout, 101) 'internal a ', 'artifical a', 'external a =', 'total a '
        write (stdout, 100) indvxdt(1, moni_particle), ardvxdt(1, moni_particle), &
            exdvxdt(1, moni_particle), dvx(1, moni_particle)
    end if

102 format(/a, i0)
101 format(4(a12,:,2x))
100 format(4(es12.5,:,2x))

end subroutine single_step
