program tester

    use test_utils, only: test_utils_to_string
    use test_param, only: test_param_keyword
    use test_input, only: test_input_unit
    use config_m, only: rk
    use test_tree_search_m, only: collect_tree_search_tests
    use parameter
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type
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
    integer stat, is
    type(testsuite_type), allocatable :: test_suites(:)
    character(*), parameter :: fmt = "('#', *(1x, a))"
    
    stat = 0
    
    test_suites = [ &
        new_testsuite("test_tree_search_m", collect_tree_search_tests) &
        ]
        
    do is = 1, size(test_suites)
        write (error_unit, fmt) "Testing:", test_suites(is)%name
        call run_testsuite(test_suites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

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
