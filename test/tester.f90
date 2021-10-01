program tester

    use test_utils, only: test_utils_to_string
    implicit none
    
    call test_utils_to_string
    print *, "All tests passed. 🎉"

end program tester