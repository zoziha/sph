!>   subroutine to calculate the smoothing kernel wij and its
!>   derivatives dwdxij.
!>     if skf = 1, cubic spline kernel by w4 - spline (monaghan 1985)
!>            = 2, gauss kernel   (gingold and monaghan 1981)
!>            = 3, quintic kernel (morris 1997)
!>
!>     r    : distance between particles i and j                     [in]
!>     dx   : x-, y- and z-distance between i and j                  [in]
!>     hsml : smoothing length                                       [in]
!>     w    : kernel for all interaction pairs                      [out]
!>     dwdx : derivative of kernel with respect to x, y and z       [out]

subroutine kernel(r, dx, hsml, w, dwdx)

    use sph_kind, only: rk
    use parameter
    implicit none

    real(rk) :: r, dx(dim), hsml, w, dwdx(dim)
    integer  :: i, j, d
    real(rk) :: q, dw, factor

    q = r/hsml
    w = 0._rk
    do d = 1, dim
        dwdx(d) = 0._rk
    end do

    if (skf == 1) then

        if (dim == 1) then
            factor = 1._rk/hsml
        else if (dim == 2) then
            factor = 15._rk/(7._rk*pi*hsml*hsml)
        else if (dim == 3) then
            factor = 3._rk/(2._rk*pi*hsml*hsml*hsml)
        else
            print *, ' >>> error <<< : wrong dimension: dim =', dim
            stop
        end if
        if (q >= 0 .and. q <= 1._rk) then
            w = factor*(2._rk/3._rk-q*q + q**3._rk/2._rk)
            do d = 1, dim
                dwdx(d) = factor*(-2._rk+3._rk/2._rk*q)/hsml**2*dx(d)
            end do
        else if (q > 1._rk .and. q <= 2) then
            w = factor*1._rk/6._rk*(2._rk-q)**3
            do d = 1, dim
                dwdx(d) = -factor*1._rk/6._rk*3.*(2._rk-q)**2/hsml*(dx(d)/r)
            end do
        else
            w = 0._rk
            do d = 1, dim
                dwdx(d) = 0._rk
            end do
        end if

    else if (skf == 2) then

        factor = 1._rk/(hsml**dim*pi**(dim/2._rk))
        if (q >= 0 .and. q <= 3) then
            w = factor*exp(-q*q)
            do d = 1, dim
                dwdx(d) = w*(-2._rk*dx(d)/hsml/hsml)
            end do
        else
            w = 0._rk
            do d = 1, dim
                dwdx(d) = 0._rk
            end do
        end if

    else if (skf == 3) then

        if (dim == 1) then
            factor = 1._rk/(120._rk*hsml)
        else if (dim == 2) then
            factor = 7._rk/(478._rk*pi*hsml*hsml)
        else if (dim == 3) then
            factor = 1._rk/(120._rk*pi*hsml*hsml*hsml)
        else
            print *, ' >>> error <<< : wrong dimension: dim =', dim
            stop
        end if
        if (q >= 0 .and. q <= 1) then
            w = factor*((3 - q)**5 - 6*(2 - q)**5 + 15*(1 - q)**5)
            do d = 1, dim
                dwdx(d) = factor*((-120 + 120*q - 50*q**2)/hsml**2*dx(d))
            end do
        else if (q > 1 .and. q <= 2) then
            w = factor*((3 - q)**5 - 6*(2 - q)**5)
            do d = 1, dim
                dwdx(d) = factor*(-5*(3 - q)**4 + 30*(2 - q)**4)/hsml*(dx(d)/r)
            end do
        else if (q > 2 .and. q <= 3) then
            w = factor*(3 - q)**5
            do d = 1, dim
                dwdx(d) = factor*(-5*(3 - q)**4)/hsml*(dx(d)/r)
            end do
        else
            w = 0._rk
            do d = 1, dim
                dwdx(d) = 0._rk
            end do
        end if

    end if

end subroutine kernel
