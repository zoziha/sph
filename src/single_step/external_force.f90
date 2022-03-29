module external_force_m
    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: ext_force
contains
    !> 计算外力的子程序，例如重力。
    !> 同时在此程序中也将边界虚粒子施加的相互作用力作为外力来计算。
    !> Subroutine to calculate the external forces, e.g. gravitational forces.
    !>  the forces from the interactions with boundary virtual particles
    !>  are also calculated here as external forces.
    subroutine ext_force(ntotal, mass, x, niac, pair_i, pair_j, itype, hsml, dvxdt)

        !> 在模拟中所使用的粒子总数
        !> number of particles in simulation
        integer, intent(in) :: ntotal
        !> 粒子的质量
        !> mass of particles
        real(rk), intent(in) :: mass(:)
        !> 粒子的坐标
        !> coordinates of particles
        real(rk), intent(in) :: x(:, :)
        !> 各粒子的粒子数
        !> number of particles in each particle
        integer, intent(in) :: niac
        !> 各粒子对的第一个粒子
        !> first particle of each particle
        integer, intent(in) :: pair_i(:)
        !> 各粒子对的第二个粒子
        !> second particle of each particle
        integer, intent(in) :: pair_j(:)
        !> 各粒子的类型
        !> type of particles
        integer, intent(in) :: itype(:)
        !> 各粒子的光滑长度
        !> smoothing length of particles
        real(rk), intent(in) :: hsml(:)
        !> 各粒子的加速度
        !> acceleration of particles
        real(rk), intent(out) :: dvxdt(:, :)

        integer :: i, j, k, d
        real(rk) :: dx(dim), rr, f, rr0, dd, p1, p2

        do i = 1, ntotal
            do d = 1, dim
                dvxdt(d, i) = 0._rk
            end do
        end do

        !     consider self-gravity or not ?

        if (self_gravity) then
            do i = 1, ntotal
                dvxdt(dim, i) = -9.8_rk
            end do
        end if

        !     boundary particle force and penalty anti-penetration force.
        rr0 = 1.25e-5_rk
        dd = 1.e-2_rk
        p1 = 12
        p2 = 4

        do k = 1, niac
            i = pair_i(k)
            j = pair_j(k)
            if (itype(i) > 0 .and. itype(j) < 0) then
                rr = 0._rk
                do d = 1, dim
                    dx(d) = x(d, i) - x(d, j)
                    rr = rr + dx(d)*dx(d)
                end do
                rr = sqrt(rr)
                if (rr < rr0) then
                    f = ((rr0/rr)**p1 - (rr0/rr)**p2)/rr**2
                    do d = 1, dim
                        dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
                    end do
                end if
            end if
        end do

    end subroutine ext_force
end module external_force_m
