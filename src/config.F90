!> SPH 精度
module config_m

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, output_unit, input_unit, error_unit
    implicit none
    private
    
    public :: rk, stdout, stdin, stderr, sp, dp
    
    !> sph预设浮点精度
    !> Default real precision for SPH
    integer, parameter :: rk = dp
    
    !> 标准输入unit
    integer, parameter :: stdin = input_unit
    !> 标准输出unit
    integer, parameter :: stdout = output_unit
    !> 标准错误unit
    integer, parameter :: stderr = error_unit

end module config_m
