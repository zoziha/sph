!> 打印当前时间和日期
!> TIME_PRINT: Print out the current date and time.
subroutine time_print()

    use sph_kind, only: rk
    use, intrinsic :: iso_fortran_env, only: stdout => output_unit
    implicit none

    ! . local scalars.
    character(len=8) :: datstr
    character(len=10) :: timstr

    ! . Get the current date and time.
    call date_and_time(datstr, timstr)

    ! . Write out the date and time.
    write (stdout, "(/A)") "                  Date = "//datstr(7:8)//"/"// &
        datstr(5:6)//"/"// &
        datstr(1:4)
    write (stdout, "(A)") "                  Time = "//timstr(1:2)//":"// &
        timstr(3:4)//":"// &
        timstr(5:10)
    write (stdout, *)

end subroutine time_print

