!> SPH 精度
module config_m

    use, intrinsic :: iso_fortran_env, only: real32, real64, output_unit, input_unit
    implicit none
    private
    
    public :: rk, stdout, stdin
    
    !> sph预设浮点精度
    !> Default real precision for SPH
    integer, parameter :: rk = real64
    
    !> 标准输入unit
    integer, parameter :: stdin = input_unit
    !> 标准输出unit
    integer, parameter :: stdout = output_unit

end module config_m
