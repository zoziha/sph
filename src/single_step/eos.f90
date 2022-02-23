!> 应用状态方程通过密度和能量计算压力的子程序
!> gamma law eos: subroutine to calculate the pressure and sound
subroutine p_gas(rho, u, p, c)

    use sph_kinds, only: rk
    implicit none
    !> 密度
    !> Density
    real(rk), intent(in) :: rho
    !> 内能
    !> Internal energy
    real(rk), intent(in) :: u
    !> 压力
    !> Pressure
    real(rk), intent(out) :: p
    !> 声速
    !> Sound velocity
    real(rk), intent(out) :: c

    real(rk) :: gamma

    !      for air (idea gas)
    !      see equ.(3.82)

    gamma = 1.4_rk
    p = (gamma - 1)*rho*u
    c = sqrt((gamma - 1)*u)

end subroutine p_gas

!> 适用人工压缩性的人工状态方程
!>   artificial equation of state for the artificial compressibility.
subroutine p_art_water(rho, p, c)

    use sph_kinds, only: rk
    implicit none
    !> 密度
    !> Density
    real(rk), intent(in) :: rho
    !> 压力
    !> Pressure
    real(rk), intent(out) :: p
    !> 声速
    !> Sound velocity
    real(rk), intent(out) :: c

    real(rk) :: gamma, rho0

    !     artificial eos, form 1 (monaghan, 1994)
    !     see equ.(4.88)
    !      gamma=7.
    !      rho0=1000.
    !      b = 1.013e5
    !      p = b*((rho/rho0)**gamma-1)
    !      c = 1480.

    !     artificial eos, form 2 (morris, 1997)
    !     see equ.(4.89)
    c = 0.01_rk
    p = c**2*rho

end subroutine p_art_water
