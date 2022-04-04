!@todo: 添加测试
module test_kernel_m

    use testdrive, only: new_unittest, unittest_type, error_type, check
    use kernel_m
    use config_m, only: rk, skf
    use parameter
    use easy_math_m, only: is_close
    implicit none
    private

    public :: collect_kernel

contains

    subroutine collect_kernel(test_suite)
        type(unittest_type), allocatable, intent(out) :: test_suite(:)

        test_suite = [ &
                     new_unittest("subroutine: kernel", test_kernel) &
                     ]

    end subroutine collect_kernel

    !@note: 也许可以把核函数类型分开测试
    subroutine test_kernel(error)
        type(error_type), allocatable, intent(out) :: error
        real(rk) :: w, dwdx(2)

        skf = 1
        call kernel(1.0_rk, [0.0_rk, 1.0_rk], 1.0_rk, w, dwdx)
        call check(error, w, 0.11368210220849664_rk, "skf1:w")
        if (allocated(error)) return
        call check(error, all(is_close(dwdx, [-0.0000000000000000_rk, &
                                              -0.34104630662549001_rk])), 'skf1:dwdx')
        if (allocated(error)) return
        
        skf = 2
        call kernel(1.0_rk, [0.0_rk, 1.0_rk], 1.0_rk, w, dwdx)
        call check(error, w, 0.11709966304863834_rk, "skf2:w")
        if (allocated(error)) return
        call check(error, all(is_close(dwdx, [-0.0000000000000000_rk, &
                                              -0.23419932609727667_rk])), 'skf2:dwdx')
        if (allocated(error)) return
        
        skf = 3
        call kernel(1.0_rk, [0.0_rk, 1.0_rk], 1.0_rk, w, dwdx)
        call check(error, w, 0.12119748804487428_rk, "skf3:w")
        if (allocated(error)) return
        call check(error, all(is_close(dwdx, [-0.0000000000000000_rk, &
                                              -0.23307209239398899_rk])), 'skf3:dwdx')

    end subroutine test_kernel

end module test_kernel_m
