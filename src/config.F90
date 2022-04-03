!> SPH 精度
module config_m

    use, intrinsic :: iso_fortran_env, only: &
        sp => real32, dp => real64, stdout => output_unit, &
        stdin => input_unit, stderr => error_unit
    implicit none
    public

    !> sph预设浮点精度, @todo: not compatible with single precision at this time
    !> Default real precision for SPH
    integer, parameter :: rk = dp

    !> 粒子搜索时间
    real, save :: tinsert = 0.0, tsearch = 0.0

    !> 输入、输出文件夹
    character(:), allocatable :: in_path !! 输入
    character(:), allocatable :: out_path !! 输出

    character(:), allocatable :: nick !! 工程名
    real(rk) :: dt !! 时间步长
    integer :: skf !! 光滑核函数的指示变量
    integer :: nnps !! 最近相邻粒子搜索算法指示变量

    !> 输出到屏幕、磁盘的时间步间隔
    integer :: print_step, save_step

end module config_m
