!>   gamma law eos: subroutine to calculate the pressure and sound
!>
!>     rho    : density                                              [in]
!>     u      : internal energy                                      [in]
!>     p      : pressure                                            [out]
!>     c      : sound velocity                                      [out]

subroutine p_gas(rho, u, p, c)

    use sph_kind, only: rk
    implicit none
    real(rk) :: rho, u, p, c
    real(rk) :: gamma

    !      for air (idea gas)
    !      see equ.(3.82)

    gamma = 1.4_rk
    p = (gamma - 1)*rho*u
    c = sqrt((gamma - 1)*u)

end subroutine p_gas

!>   artificial equation of state for the artificial compressibility.
!>
!>     rho    : density                                              [in]
!>     u      : internal energy                                      [in]
!>     p      : pressure                                            [out]
!>     c      : sound velocity                                      [out]
!>     equation of state for artificial compressibility

subroutine p_art_water(rho, p, c)

    use sph_kind, only: rk
    implicit none
    real(rk) :: rho, u, p, c
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
