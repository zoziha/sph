module sph_kind

    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private

    public :: rk

    !> sph预设浮点精度
    !> Default real precision for SPH
    integer, parameter :: rk = real64

end module sph_kind
