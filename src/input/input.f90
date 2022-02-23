!> 装载或产生初始数据的子程序。
!> subroutine for loading or generating initial particle information
subroutine input(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    use sph_kinds, only: rk
    use parameter
    implicit none

    !> 粒子的位置
    !> coordinates of particles
    real(rk), intent(out) :: x(dim, maxn)
    !> 粒子的速度
    !> velocities of particles
    real(rk), intent(out) :: vx(dim, maxn)
    !> 粒子的质量
    !> mass of particles
    real(rk), intent(out) :: mass(maxn)
    !> 粒子的密度
    !> dnesities of particles
    real(rk), intent(out) :: rho(maxn)
    !> 粒子的压力
    !> pressure  of particles
    real(rk), intent(out) :: p(maxn)
    !> 粒子的内部能量
    !> internal energy of particles
    real(rk), intent(out) :: u(maxn)
    !> 粒子的类型
    !> types of particles
    integer, intent(out) :: itype(maxn)
    !> 粒子的光滑长度
    !> smoothing lengths of particles
    real(rk), intent(out) :: hsml(maxn)
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(out) :: ntotal

    integer :: i, d, im

    !> load initial particle information from external disk file

    if (config_input) then

        open (1, file='./data/f_xv.dat')
        open (2, file='./data/f_state.dat')
        open (3, file='./data/f_other.dat')

        write (*, *) '  **************************************************'
        write (*, *) '      loading initial particle configuration...   '
        read (1, *) ntotal
        write (*, '(1x,a,i0)') '      total number of particles: ', ntotal
        write (*, *) '  **************************************************'
        do i = 1, ntotal
            read (1, *) im, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            read (2, *) im, mass(i), rho(i), p(i), u(i)
            read (3, *) im, itype(i), hsml(i)
        end do

    else

        open (1, file='./data/ini_xv.dat')
        open (2, file='./data/ini_state.dat')
        open (3, file='./data/ini_other.dat')

        if (shocktube) call shock_tube(x, vx, mass, rho, p, u, itype, hsml, ntotal)

        if (shearcavity) call shear_cavity(x, vx, mass, rho, p, u, itype, hsml, ntotal)
        do i = 1, ntotal
            write (1, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            write (2, 1002) i, mass(i), rho(i), p(i), u(i)
            write (3, 1003) i, itype(i), hsml(i)
        end do
        write (*, *) '  **************************************************'
        write (*, *) '      initial particle configuration generated   '
        write (*, '(1x,a,i0)') '      total number of particles: ', ntotal
        write (*, *) '  **************************************************'

    end if

    close (1)
    close (2)
    close (3)

1001 format(1x, i5, 6(2x, e14.8))
1002 format(1x, i5, 7(2x, e14.8))
1003 format(1x, i5, 2x, i2, 2x, e14.8)

end subroutine input

!> 一维振荡管的初始数据。
!> this subroutine is used to generate initial data for the
!> 1 d noh shock tube problem
subroutine shock_tube(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    use sph_kinds, only: rk
    use parameter
    implicit none

    !> 粒子的位置
    !> coordinates of particles
    real(rk), intent(out) :: x(dim, maxn)
    !> 粒子的速度
    !> velocities of particles
    real(rk), intent(out) :: vx(dim, maxn)
    !> 粒子的质量
    !> mass of particles
    real(rk), intent(out) :: mass(maxn)
    !> 粒子的密度
    !> dnesities of particles
    real(rk), intent(out) :: rho(maxn)
    !> 粒子的压力
    !> pressure  of particles
    real(rk), intent(out) :: p(maxn)
    !> 粒子的内部能量
    !> internal energy of particles
    real(rk), intent(out) :: u(maxn)
    !> 粒子的类型(1: ideal gas; 2: water)
    !> types of particles
    integer, intent(out) :: itype(maxn)
    !> 粒子的光滑长度
    !> smoothing lengths of particles
    real(rk), intent(out) :: hsml(maxn)
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(out) :: ntotal

    integer :: i, d
    real(rk) :: space_x

    ntotal = 400
    space_x = 0.6_rk/80._rk

    do i = 1, ntotal
        mass(i) = 0.75_rk/400._rk
        hsml(i) = 0.015_rk
        itype(i) = 1
        do d = 1, dim
            x(d, i) = 0._rk
            vx(d, i) = 0._rk
        end do
    end do

    do i = 1, 320
        x(1, i) = -0.6_rk + space_x/4._rk*(i - 1)
    end do

    do i = 320 + 1, ntotal
        x(1, i) = 0._rk + space_x*(i - 320)
    end do

    do i = 1, ntotal
        if (x(1, i) <= 1.e-8_rk) then
            u(i) = 2.5_rk
            rho(i) = 1._rk
            p(i) = 1._rk
        end if
        if (x(1, i) > 1.e-8_rk) then
            u(i) = 1.795_rk
            rho(i) = 0.25_rk
            p(i) = 0.1795_rk
        end if
    end do

end subroutine shock_tube

!> 二维剪切腔的初始数据。
!> this subroutine is used to generate initial data for the
!> 2 d shear driven cavity probem with re = 1
subroutine shear_cavity(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    use sph_kinds, only: rk
    use parameter
    implicit none

    !> 粒子的位置
    !> coordinates of particles
    real(rk), intent(out) :: x(dim, maxn)
    !> 粒子的速度
    !> velocities of particles
    real(rk), intent(out) :: vx(dim, maxn)
    !> 粒子的质量
    !> mass of particles
    real(rk), intent(out) :: mass(maxn)
    !> 粒子的密度
    !> dnesities of particles
    real(rk), intent(out) :: rho(maxn)
    !> 粒子的压力
    !> pressure  of particles
    real(rk), intent(out) :: p(maxn)
    !> 粒子的内部能量
    !> internal energy of particles
    real(rk), intent(out) :: u(maxn)
    !> 粒子的类型(1: ideal gas; 2: water)
    !> types of particles
    integer, intent(out) :: itype(maxn)
    !> 粒子的光滑长度
    !> smoothing lengths of particles
    real(rk), intent(out) :: hsml(maxn)
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(out) :: ntotal

    integer :: i, j, d, m, n, mp, np, k
    real(rk) :: xl, yl, dx, dy

    !> giving mass and smoothing length as well as other data.

    m = 41
    n = 41
    mp = m - 1
    np = n - 1
    ntotal = mp*np
    xl = 1.e-3_rk
    yl = 1.e-3_rk
    dx = xl/mp
    dy = yl/np

    do i = 1, mp
        do j = 1, np
            k = j + (i - 1)*np
            x(1, k) = (i - 1)*dx + dx/2._rk
            x(2, k) = (j - 1)*dy + dy/2._rk
        end do
    end do

    do i = 1, mp*np
        vx(1, i) = 0._rk
        vx(2, i) = 0._rk
        rho(i) = 1000._rk
        mass(i) = dx*dy*rho(i)
        p(i) = 0._rk
        u(i) = 357.1_rk
        itype(i) = 2
        hsml(i) = dx
    end do

end subroutine shear_cavity
