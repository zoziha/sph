! 支持核函数: 
! 1. 三次样条核函数 (monaghan 1985);
! 2. 高斯核函数 (gingold and monaghan 1981);
! 3. 五次核函数 (morris 1997)
module kernel_m

    use config_m, only: rk, skf
    use parameter, only: pi, dim
    implicit none
    private

    public :: kernel

contains

    !> 计算光滑函数 Wij 及其导数 dWdxij 的子例程
    pure subroutine kernel(r, dx, hsml, w, dwdx)
        real(rk), intent(in) :: r, dx(dim), hsml
        real(rk), intent(out) :: w, dwdx(dim)
        real(rk) :: q, factor

        q = r/hsml

        if (skf == 1) then

            if (dim == 1) then
                factor = 1.0_rk/hsml
            elseif (dim == 2) then
                factor = 15._rk/(7.0_rk*pi*hsml*hsml)
            elseif (dim == 3) then
                factor = 3.0_rk/(2.0_rk*pi*hsml**3)
            end if

            if (q >= 0 .and. q <= 1) then
                w = factor*(2.0_rk/3.0_rk - q*q + q**3/2.0_rk)
                dwdx = factor*(-2.0_rk + 3.0_rk/2.0_rk*q)/hsml**2*dx
            elseif (q > 1 .and. q <= 2) then
                w = factor/6.0_rk*(2.0_rk - q)**3
                dwdx = -factor/6.0_rk*3.0_rk*(2.0_rk - q)**2/hsml*(dx/r)
            else
                w = 0.0_rk
                dwdx = 0.0_rk
            end if

        elseif (skf == 2) then

            factor = 1.0_rk/(hsml**dim*pi**(dim/2.0_rk))
            if (q >= 0 .and. q <= 3) then
                w = factor*exp(-q*q)
                dwdx = w*(-2.0_rk*dx/(hsml*hsml))
            else
                w = 0.0_rk
                dwdx = 0.0_rk
            end if

        elseif (skf == 3) then

            if (dim == 1) then
                factor = 1.0_rk/(120.0_rk*hsml)
            elseif (dim == 2) then
                factor = 7.0_rk/(478.0_rk*pi*hsml*hsml)
            elseif (dim == 3) then
                factor = 1.0_rk/(120.0_rk*pi*hsml**3)
            end if

            if (q >= 0 .and. q <= 1) then
                w = factor*((3 - q)**5 - 6*(2 - q)**5 + 15*(1 - q)**5)
                dwdx = factor*((-120 + 120*q - 50*q**2)/hsml**2*dx)
            elseif (q > 1 .and. q <= 2) then
                w = factor*((3 - q)**5 - 6*(2 - q)**5)
                dwdx = factor*(-5*(3 - q)**4 + 30*(2 - q)**4)/hsml*(dx/r)
            elseif (q > 2 .and. q <= 3) then
                w = factor*(3 - q)**5
                dwdx = factor*(-5*(3 - q)**4)/hsml*(dx/r)
            else
                w = 0.0_rk
                dwdx = 0.0_rk
            end if

        end if

    end subroutine kernel

end module kernel_m
