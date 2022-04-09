!> 状态方程
!>
!> 1. 气体;
!> 2. 淡水
module eos_m

    use config_m, only: rk, eos_form
    implicit none
    private

    public :: p_art_water, p_gas

contains

    !> 应用状态方程通过密度和能量计算压力
    pure subroutine p_gas(rho, u, p, c)
        real(rk), intent(in) :: rho !! 密度
        real(rk), intent(in) :: u   !! 内能
        real(rk), intent(out) :: p  !! 压力
        real(rk), intent(out) :: c  !! 声速

        real(rk), parameter :: gamma = 1.4_rk

        ! for air (idea gas)
        ! see equ.(3.82)
        p = (gamma - 1)*rho*u
        c = sqrt((gamma - 1)*u)

    end subroutine p_gas

    !> 适用人工压缩性的人工状态方程
    pure subroutine p_art_water(rho, p, c)
        real(rk), intent(in) :: rho !! 密度
        real(rk), intent(out) :: p  !! 压力
        real(rk), intent(out) :: c  !! 声速

        real(rk) :: rho0, b
        integer, parameter :: gamma = 7

        select case (eos_form)
        case (1)
            ! artificial eos, form 1 (monaghan, 1994)
            ! see equ.(4.88)
            rho0 = 1000.0_rk
            b = 1.013e5_rk          ! B = rho0 * c^2 / gamma
            p = b*((rho/rho0)**gamma - 1)
            c = 1480.0_rk
        case (2)
            ! artificial eos, form 2 (morris, 1997)
            ! see equ.(4.89)
            c = 0.01_rk
            p = c**2*rho
        end select
    end subroutine p_art_water

end module eos_m
