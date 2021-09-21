!>   subroutine to define the fluid particle viscosity
!>
!>     ntotal  : number of particles                                 [in]
!>     itype    : type of particle                                   [in]
!>     x       : coordinates of all particles                        [in]
!>     rho     : density                                             [in]
!>     eta     : dynamic viscosity                                  [out]

subroutine viscosity(ntotal, itype, x, rho, eta)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: ntotal, i, itype(maxn)
    real(rk) :: x(dim, maxn), rho(maxn), eta(maxn)

    do i = 1, ntotal
        if (abs(itype(i)) == 1) then
            eta(i) = 0._rk
        else if (abs(itype(i)) == 2) then
            eta(i) = 1.0e-3_rk
        end if
    end do

end subroutine viscosity
