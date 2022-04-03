!@todo: 添加测试
module test_macro_m

    use testdrive, only: new_unittest, unittest_type, error_type, check
    use macro_m
    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: collect_macro

contains

    subroutine collect_macro(test_suite)
        type(unittest_type), allocatable, intent(out) :: test_suite(:)

        test_suite = [ &
                     new_unittest("subroutine: alloc_macro_memory", test_alloc_macro_memory) &
                     ]

    end subroutine collect_macro

    subroutine test_alloc_macro_memory(error)
        type(error_type), allocatable, intent(out) :: error
        
        call alloc_macro_memory(1, 3, x, vx, mass, rho, p, u, c, s, e, hsml, itype)
        call check(error, size(x, 1), 1, "alloc_macro_memory: x, 1")
        if (allocated(error)) return
        call check(error, size(x, 2), 3, "alloc_macro_memory: x, 2")

    end subroutine test_alloc_macro_memory

end module test_macro_m
