!> SPH 精度
module config_m

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, output_unit, input_unit, error_unit
    implicit none
    public
    
    !> sph预设浮点精度
    !> Default real precision for SPH
    integer, parameter :: rk = dp
    
    !> 标准输入unit
    integer, parameter :: stdin = input_unit
    !> 标准输出unit
    integer, parameter :: stdout = output_unit
    !> 标准错误unit
    integer, parameter :: stderr = error_unit
    
    !> 粒子搜索时间
    real, save :: tinsert = 0.0, tsearch = 0.0
    
    !> 输入、输出文件夹
    character(:), allocatable :: in_path !! 输入
    character(:), allocatable :: out_path !! 输出
    
    character(:), allocatable :: nick !! 工程名
    real(rk) :: dt !! 时间步长
    integer :: skf !! 光滑核函数的指示变量
    integer :: nnps !! 最近相邻粒子搜索算法指示变量

end module config_m
