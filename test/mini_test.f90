module mini_test

    implicit none

contains

    subroutine check(condition, msg)
        logical,intent(in)           :: condition
        character(len=*), intent(in) :: msg

        if (.not.condition) error stop msg

    end subroutine check

end module mini_test