!> 读取命令行
module arg_m

    use argv_module, only: argv
    implicit none
    private

    public :: get_path

contains

    !> 获取第一个命令行参数为路径名
    subroutine get_path(path)
        character(:), allocatable, intent(out) :: path !! 路径名
        integer :: nargs
        nargs = command_argument_count()
        if (nargs < 1) then
            path = '.'
        else
            path = argv(1)
        end if
    end subroutine get_path

end module arg_m
