!> 通用工具
!> General Tools
module utils

    use sph_kinds, only: rk
    use parameter
    implicit none
    private

    public :: to_string, tic, toc, get_distance

    integer, save :: time_save

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
        call system_clock(time_save)
    end subroutine tic

    subroutine toc()
        integer :: time_now
        call system_clock(time_now)
        write (*, "(A, i0, A)") "Elapsed cpu time (seconds) = ", &
            (time_now - time_save)/1000, "s"
    end subroutine toc

    !> 获取两点之间的距离
    subroutine get_distance(x, y, d, r)
        real(rk), intent(in), dimension(dim) :: x, y
        real(rk), intent(out), dimension(dim) :: d
        real(rk), intent(out) :: r
        integer i

        d(1) = x(1) - y(1)
        r = d(1)*d(1)

        do i = 2, dim
            d(i) = x(i) - y(i)
            r = r + d(i)*d(i)
        end do

    end subroutine get_distance

end module utils
