!> 通用工具
!> General Tools
module utils

    use sph_kind, only: rk
    implicit none
    private
    
    public :: to_string, tic, toc
    
    real(rk), save :: time_save

contains
    
    ! SPDX-Identifier: MIT
    !> 将整数转化为字符串。借鉴了Fortran标准库的to_string。
    !> Change integer to string.
    pure function to_string(x) result(string)

        !> 需要转换为字符串的整数
        !> Integer to be converted to string
        integer, intent(in) :: x
        !> 转换后的字符串
        !> String after conversion
        character(:), allocatable :: string

        integer, parameter :: buffer_len = range(x) + 2
        character(buffer_len) :: buffer
        integer :: pos
        integer :: n
        character(1), parameter :: number(0:9) = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]

        if (x == 0) then
            string = number(0)
            return
        end if

        n = abs(x)
        buffer = ""

        pos = buffer_len + 1
        do while (n > 0)
            pos = pos - 1
            buffer(pos:pos) = number(mod(n, 10))
            n = n/10
        end do
        if (x < 0) then
            pos = pos - 1
            buffer(pos:pos) = "-"
        end if

        string = buffer(pos:)

    end function to_string
    
    ! 非并行方案
    subroutine tic()
        call cpu_time(time_save)
    end subroutine tic

    subroutine toc()
        real(rk) :: time_now
        call cpu_time(time_now)
        write (*, "(A, F0.1)") "elapsed cpu time (seconds) = ", time_now - time_save
    end subroutine toc

end module utils
