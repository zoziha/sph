module mini_test

    use config_m, only: rk
    implicit none
    private

    public :: check

    !> 断言
    interface check
        procedure check_logial
        procedure check_real
    end interface check

contains

    !> 断言
    subroutine check_logial(condition, msg)
        logical, intent(in) :: condition     !! 断言条件
        character(len=*), intent(in) :: msg  !! 错误提示

        ! 需要Fortran 2013标准以上，以支持error stop
        if (.not. condition) error stop msg

    end subroutine check_logial

    !> 断言
    subroutine check_real(r1, r2, msg)
        real(rk), intent(in) :: r1  !! 第一个数
        real(rk), intent(in) :: r2  !! 第二个数
        character(len=*), intent(in) :: msg  !! 错误提示

        if (abs(r1 - r2) >= 1.0e-3_rk) error stop msg

    end subroutine check_real

end module mini_test
