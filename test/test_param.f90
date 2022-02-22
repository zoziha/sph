module test_param

    use mini_test
    use parameter
    implicit none

contains

    !> 测试parameter模块中的关键变量
    subroutine test_param_keyword()

        call check(abs(pi - 3.14159265358979) <= 1.0e-7, msg="pi 值不正确")

    end subroutine test_param_keyword

end module test_param
