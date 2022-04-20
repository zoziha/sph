!@todo: 添加测试
module test_output_m

    use testdrive, only: new_unittest, unittest_type, error_type, check
    use output_m
    use config_m, only: rk, skf, print_step
    use parameter
    implicit none
    private

    public :: collect_output

contains

    subroutine collect_output(test_suite)
        type(unittest_type), allocatable, intent(out) :: test_suite(:)

        test_suite = [ &
                     new_unittest("sub: output", test_set_statistics_print) &
                     ]

    end subroutine collect_output

    subroutine test_set_statistics_print(error)
        type(error_type), allocatable, intent(out) :: error
        
        skf = 1
        print_step = 100
        call set_statistics_print(100, 4, 2, [1, 1, 0, 0])

    end subroutine test_set_statistics_print

end module test_output_m
