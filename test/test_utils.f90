module test_utils

    use utils, only: to_string
    use mini_test, only: check

contains

    subroutine test_utils_to_string

        call check(to_string(12)=="12", msg='to_string(12)=="12" failed.')

    end subroutine test_utils_to_string

end module test_utils