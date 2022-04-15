!@todo: 添加测试
module test_toml_info_m

    use testdrive, only: new_unittest, unittest_type, error_type, check
    use toml_info_m, only: parse_toml_info
    use config_m
    implicit none
    private

    public :: collect_toml_info

contains

    subroutine collect_toml_info(test_suite)
        type(unittest_type), allocatable, intent(out) :: test_suite(:)

        test_suite = [ &
                     new_unittest("subroutine: parse_toml_info", test_parse_toml_info) &
                     ]

    end subroutine collect_toml_info

    subroutine test_parse_toml_info(error)
        type(error_type), allocatable, intent(out) :: error
        call parse_toml_info()
        call check(error, in_path, './data/demo/shearcavity')
        if (allocated(error)) return
        call check(error, out_path, './data/output')
        if (allocated(error)) return
        call check(error, nick, "2 dimensional shear cavity flow")
        if (allocated(error)) return
        call check(error, dt, 5.0e-5_rk)
        if (allocated(error)) return
        call check(error, skf, 1)
        if (allocated(error)) return
        call check(error, nnps, 2)
        if (allocated(error)) return
        call check(error, print_step, 100)
        if (allocated(error)) return
        call check(error, save_step, 500)
        if (allocated(error)) return
        call check(error, maxn, 3000)
        if (allocated(error)) return
        call check(error, max_interaction, maxn*20)
        if (allocated(error)) return

    end subroutine test_parse_toml_info

end module test_toml_info_m
