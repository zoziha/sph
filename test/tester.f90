program tester

    use, intrinsic :: iso_fortran_env, only: error_unit
    use config_m, only: rk
    use parameter
    use stdlib_logger, only: stdlog => global_logger
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type
    use test_eos_m, only: collect_eos
    use test_macro_m, only: collect_macro
    use test_kernel_m, only: collect_kernel
    use test_output_m, only: collect_output
    use test_toml_info_m, only: collect_toml_info
    use test_tree_search_m, only: collect_tree_search_tests
    implicit none
    integer stat, is
    type(testsuite_type), allocatable :: test_suites(:)
    character(*), parameter :: fmt = "('#', *(1x, a))"

    call stdlog%add_log_file('test/.stdlog.log')    !! 输出日志文件
    call stdlog%log_information('Starting testing')
    stat = 0
    !@todo: test art_heat/art_visc
    test_suites = [ &
                  new_testsuite("mod: test_toml_info_m", collect_toml_info), &
                  new_testsuite("mod: test_macro_m", collect_macro), &
                  new_testsuite("mod: test_kernel_m", collect_kernel), &
                  new_testsuite("mod: test_eos_m", collect_eos), &
                  new_testsuite("mod: test_tree_search_m", collect_tree_search_tests) &
                  ]

    do is = 1, size(test_suites)
        write (error_unit, fmt) "Testing:", test_suites(is)%name
        call run_testsuite(test_suites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

    print '(/a)', "# ******** All tests passed :)"
    call stdlog%log_information('Testing finished')

end program tester
