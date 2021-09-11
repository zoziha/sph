   subroutine time_elapsed(s)

!===============================================================================
!   The standard Fortran 90 routine RTC is used to calculate the elapsed CPU
!===============================================================================

   !  use dfport
    implicit none

    integer, parameter :: output = 6
    real(8) :: s

   ! s = rtc()
    s = 1._8

   end subroutine time_elapsed