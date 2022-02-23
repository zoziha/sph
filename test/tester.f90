program tester

    use test_utils, only: test_utils_to_string
    use test_param, only: test_param_keyword
    use test_input, only: test_input_unit
    use sph_kinds, only: rk
    use parameter
    implicit none

    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer :: ntotal
    !> 粒子的类型(1: ideal gas; 2: water)
    !> types of particles
    integer :: itype(maxn)
    integer :: maxtimestep
    integer :: d, m, i, yesorno
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), &
                u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), dt

    ! 辅助单元
    call test_utils_to_string()
    call test_param_keyword()

    ! 业务单元
    if (shocktube) dt = 0.005_rk
    if (shearcavity) dt = 5.0e-5_rk
    call input(x, vx, mass, rho, p, u, itype, hsml, ntotal)
    call test_input_unit(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    print *, "All tests passed. :)"

end program tester
