!> 输出模块(TODO: 待更名)
module output_m

    use utils, only: to_string
    implicit none

contains

    !> 输出每个保存时间步的求解信息（拓展）
    !> subroutine for saving particle information to external disk file
    subroutine output_all(x, vx, mass, rho, p, u, c, itype, hsml, ntotal, n)

        use sph_kind, only: rk
        use parameter
        implicit none

        !> 粒子的位置
        !> coordinates of particles
        real(rk), intent(in) :: x(dim, maxn)
        !> 粒子的速度
        !> velocities of particles
        real(rk), intent(in) :: vx(dim, maxn)
        !> 粒子的质量
        !> mass of particles
        real(rk), intent(in) :: mass(maxn)
        !> 粒子的密度
        !> dnesities of particles
        real(rk), intent(in) :: rho(maxn)
        !> 粒子的压力
        !> pressure  of particles
        real(rk), intent(in) :: p(maxn)
        !> 粒子的内部能量
        !> internal energy of particles
        real(rk), intent(in) :: u(maxn)
        !> 粒子的声速
        !> sound velocity of particles
        real(rk), intent(in) :: c(maxn)
        !> 粒子的类型(1: ideal gas; 2: water; 3: TNT)
        !> types of particles
        integer, intent(in) :: itype(maxn)
        !> 粒子的光滑长度
        !> smoothing lengths of particles
        real(rk), intent(in) :: hsml(maxn)
        !> 在模拟中所使用的粒子总数
        !> number of particles in simulation
        integer, intent(in) :: ntotal
        !> 第几个时间步
        !> number of time step
        integer, intent(in) :: n

        integer :: i, d, npart
        integer :: xv_unit, state_unit, other_unit

        !> 输出粒子的位置、速度信息
        open (newunit=xv_unit, file='./data/all/f_'//to_string(n)//'xv.dat')

        !> 输出粒子的宏观信息：质量、密度、压强、内能
        open (newunit=state_unit, file='./data/all/f_'//to_string(n)//'state.dat')

        !> 输出粒子的其它信息：粒子类型、光滑长度
        open (newunit=other_unit, file='./data/all/f_'//to_string(n)//'other.dat')

        write (xv_unit, *) ntotal
        do i = 1, ntotal
            write (xv_unit, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            write (state_unit, 1002) i, mass(i), rho(i), p(i), u(i)
            write (other_unit, 1003) i, itype(i), hsml(i)
        end do

        close (xv_unit)
        close (state_unit)
        close (other_unit)

1001    format(1x, i6, 6(2x, e14.8))
1002    format(1x, i6, 7(2x, e14.8))
1003    format(1x, i6, 2x, i4, 2x, e14.8)

    end subroutine output_all

end module output_m
