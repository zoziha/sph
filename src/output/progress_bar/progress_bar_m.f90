!> 进度条
!> char(8): 退格
!> char(13): 回车
!> 借鉴、修改参考: https://github.com/StellaContrail/FortranProgressbar
!> 如有侵权，删除此模块，并向作者 (zuo.zhihua@qq.com) 发送邮件，谢谢！
module progress_bar_m

    use utils, only: to_string
    use sph_kinds, only: rk
    implicit none
    integer, private :: fdigit = 0

contains

    !> 利用system_clock获取运行时间
    !> 在所有输出处理 (Any) 完成后调用
    !> value / max:%（超出最大值）
    !> isflushed：当有计划flush时使用，主要用于在输出进度条的同时输出结果，不使用时设置为false。
    subroutine pbout(value, max, isflushed)
        integer, intent(in) :: value, max
        logical, intent(in) :: isflushed

        integer, parameter :: digit = 50
        real(rk), save :: rate = 0.0_rk
        integer, save :: time = 0
        real(rk) dr
        integer estimate
        integer remain(3), dt

        if (isflushed) then
            write (*, *)
        end if

        ! 获取日期和时间
        dr = rate   ! 将之前执行的速率暂时存储在dr中
        dt = time   ! 将之前执行的时间暂时存入dt

        ! 获取当前时间标记（毫秒）
        call system_clock(time)
        dt = time - dt

        ! 百分比更新
        rate = real(value, rk)/real(max, rk)
        dr = rate - dr

        ! 剩余时间的计算（毫秒）
        estimate = int((1.0_rk - rate)*(real(dt, rk)/dr))

        ! 转换为 hour / min / sec 表示法
        remain(1) = estimate/3600000 ! 小时
        remain(2) = (estimate/60000 - remain(1)*60) ! 分
        remain(3) = (estimate/1000 - remain(1)*3600 - remain(2)*60) ! 秒

        write (*, "(2(a, i0), a)", advance='no') "(", value, "/", max, ") ["
        write (*, '(2a)', advance='no') repeat("#", int(digit*rate)), repeat(" ", digit - int(digit*rate))

        if (isflushed) fdigit = 2*int(log10(real(max))) + 90
        write (*, '(a, f7.2, a, i3, a, 2(i2, a), a)', advance='no') "]", 100*rate, "% 剩余时间:", &
            remain(1), "h", remain(2), "m", remain(3), "s", char(13)

        if (value == max) write (*, *)

    end subroutine pbout

    !> 在第一个输出处理 (Any) 开始之前调用
    !> 仅当使用 pbout () 将 isflushed 设置为 true 时才需要。 false 时不需要。
    subroutine pbflush()
        character(len=10) FMT
        if (fdigit == 0) return

        FMT = "("//to_string(fdigit)//"x, 3a)"
        write (*, FMT, advance='no') char(13), char(8), char(13)
    end subroutine pbflush

end module progress_bar_m

!> Demo
! program main
!     use progress_bar_m
!     implicit none
!     integer i

!     do i = 1, 100
!         if (mod(i, 10) == 0) then
!             call pbflush()
!             !write (*, *) "#", i
!             write (*, *) "It is a line."
!             call pbout(i, 100, .true.)
!             call sleep(1)
!         end if
!     end do

!     print *, "done"

! end program
