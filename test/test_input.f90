module test_input

    use sph_kind, only: rk
    use parameter
    use mini_test
    implicit none

contains

    subroutine test_input_unit(x, vx, mass, rho, p, u, itype, hsml, ntotal)
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
        !> 粒子的类型
        !> types of particles
        integer, intent(in) :: itype(maxn)
        !> 粒子的光滑长度
        !> smoothing lengths of particles
        real(rk), intent(in) :: hsml(maxn)
        !> 在模拟中所使用的粒子总数
        !> number of particles in simulation
        integer, intent(in) :: ntotal

        if (shearcavity) then
            ! 目前只检查粒子的少量特征信息
            call check(ntotal == 1600, "ntotal不等于1600")
            call check(x(1, 1), 0.125e-4_rk, "x(1,1)不等于0.125e-4_rk")
            call check(x(2, 1), 0.125e-4_rk, "x(2,1)不等于0.125e-4_rk")
            call check(vx(1, 1), 0.0_rk, "vx(1,1)不等于0.0_rk")
            call check(vx(2, 1), 0.0_rk, "vx(2,1)不等于0.0_rk")

            call check(mass(1), 0.625e-6_rk, "mass(1)不等于0.625e-6_rk")
            call check(rho(1), 0.1e4_rk, "rho(1)不等于0.1e4_rk")
            call check(p(1), 0.0_rk, "p(1)不等于0.0_rk")
            call check(u(1), 0.3571e3_rk, "u(1)不等于0.3571e3_rk")

            call check(itype(1) == 2, "itype(1)不等于2")
            call check(hsml(1), 0.25e-4_rk, "hsml(1)不等于0.25e-4_rk")
        end if

    end subroutine test_input_unit

end module test_input
