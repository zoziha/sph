!> 计算光滑函数 Wij 及其导数 dWdxij 的子例程
!> subroutine to calculate the smoothing kernel wij and its
!> derivatives dwdxij.
subroutine kernel(r, dx, hsml, w, dwdx)

    use sph_kind, only: rk
    use parameter
    implicit none

    !> 粒子 i 和 j 之间的距离
    !> distance between particles i and j
    real(rk), intent(in) :: r
    !> 粒子 i 和 j 之间的坐标差
    !> x-, y- and z-distance between i and j
    real(rk), intent(in) :: dx(dim)
    !> 粒子的光滑长度
    !> smoothing length
    real(rk), intent(in) :: hsml
    !> 给定相互作用对的光滑核函数
    !> kernel for all interaction pairs
    real(rk), intent(out) :: w
    !> 核函数对 x, y, z 的导数
    !> derivative of kernel with respect to x, y and z
    real(rk), intent(out) :: dwdx(dim)

    integer :: i, j, d
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
            w = factor*(2._rk/3._rk - q*q + q**3._rk/2._rk)
            do d = 1, dim
                dwdx(d) = factor*(-2._rk + 3._rk/2._rk*q)/hsml**2*dx(d)
            end do
        else if (q > 1._rk .and. q <= 2) then
            w = factor*1._rk/6._rk*(2._rk - q)**3
            do d = 1, dim
                dwdx(d) = -factor*1._rk/6._rk*3.*(2._rk - q)**2/hsml*(dx(d)/r)
                !dwda(d) =-factor*(2._rk - q)**2/2*hsml*(dx(d)/r)
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
