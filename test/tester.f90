program tester

    use config_m, only: rk
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type
    use parameter
    use, intrinsic :: iso_fortran_env, only: error_unit
    use test_tree_search_m, only: collect_tree_search_tests
    use test_macro_m, only: collect_macro
    use test_kernel_m, only: collect_kernel
    use test_output_m, only: collect_output
    implicit none
    integer stat, is
    type(testsuite_type), allocatable :: test_suites(:)
    character(*), parameter :: fmt = "('#', *(1x, a))"
    
    stat = 0
    !@todo: test art_heat/art_visc
    test_suites = [ &
        new_testsuite("module: test_tree_search_m", collect_tree_search_tests), &
        new_testsuite("module: test_macro_m", collect_macro), &
        new_testsuite("module: test_kernel_m", collect_kernel) &
        ]
        
    do is = 1, size(test_suites)
        write (error_unit, fmt) "Testing:", test_suites(is)%name
        call run_testsuite(test_suites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

    print *, "All tests passed. :)"

end program tester
