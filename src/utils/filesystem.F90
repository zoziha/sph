module filesystem

    implicit none

contains

    !> 判断当前系统是否为 Windows 系统
    logical function is_windows()
        character(16) :: os_name
        logical, save :: is_windows_ = .false.
        logical, save :: is_first_run = .true.

        if (is_first_run) then
            call get_environment_variable("OS", os_name)
            is_windows_ = trim(os_name) == "Windows_NT"
            is_first_run = .false.
            is_windows = is_windows_
        else
            is_windows = is_windows_
        end if

    end function is_windows

    !> 创建目录
    subroutine mkdir(path)
        character(*), intent(in) :: path  !! 目录路径

        if (is_windows()) then
            call execute_command_line("md "//path)
        else
            call execute_command_line("mkdir -p "//path)
        end if

    end subroutine mkdir

    !> 判断目录是否存在
    function exists(path, is_directory)
        character(*), intent(in) :: path       !! 目录路径
        logical, intent(in), optional :: is_directory    !! 是否为目录
        logical :: exists

#if defined __INTEL_COMPILER
        if (present(is_directory)) then
            if (is_directory) then
                inquire (directory=path, exist=exists)
                return
            end if
        end if
#endif
        inquire (file=path, exist=exists)

    end function exists

end module filesystem
