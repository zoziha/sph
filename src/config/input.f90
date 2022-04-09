!> 输入模块，包含:
!>
!> 1. 虚粒子信息临时存储;
!> 2. 实、虚粒子读入程序;
!> 3. 两个个典型例子，一维冲击管，二维剪切腔
module input_m

    use config_m, only: rk, stdout, in_path, save_step, print_step
    use parameter
    use info_m, only: operator(.c.)
    implicit none
    private

    public :: input, shock_tube, shear_cavity, virt_part
    public :: saved_virt_part ! public for lua_call_m

    !> 临时存储的虚拟粒子信息
    type virt_part_info_t
        integer :: nvirt
        real(rk), allocatable :: x(:, :), vx(:, :), mass(:), rho(:), p(:), hsml(:), u(:)
        integer, allocatable :: itype(:)
    end type virt_part_info_t

    !> 临时存储的虚拟粒子信息
    type(virt_part_info_t), save :: saved_virt_part

contains

    !> 装载或产生初始数据的子程序
    subroutine input(x, vx, mass, rho, p, u, itype, hsml, ntotal)
        real(rk), intent(out) :: x(:, :)    !! 粒子位置
        real(rk), intent(out) :: vx(:, :)   !! 粒子速度
        real(rk), intent(out) :: mass(:)    !! 粒子质量
        real(rk), intent(out) :: rho(:)     !! 粒子密度
        real(rk), intent(out) :: p(:)       !! 粒子压力
        real(rk), intent(out) :: u(:)       !! 粒子内部能量
        integer, intent(out) :: itype(:)    !! 粒子类型
        real(rk), intent(out) :: hsml(:)    !! 粒子光滑长度
        integer, intent(out) :: ntotal      !! 粒子总数

        integer :: i, d, im

        ! load initial particle information from external disk file

        if (config_input) then

            !@todo: in_path
            open (1, file=in_path//'/f_xv.dat')
            open (2, file=in_path//'/f_state.dat')
            open (3, file=in_path//'/f_other.dat')

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

            open (1, file=in_path//'/ini_xv.dat')
            open (2, file=in_path//'/ini_state.dat')
            open (3, file=in_path//'/ini_other.dat')

            if (shocktube) call shock_tube(x, vx, mass, rho, p, u, itype, hsml, ntotal)

            if (shearcavity) call shear_cavity(x, vx, mass, rho, p, u, itype, hsml, ntotal)
            do i = 1, ntotal
                write (1, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
                write (2, 1002) i, mass(i), rho(i), p(i), u(i)
                write (3, 1003) i, itype(i), hsml(i)
            end do
            write (stdout, '(a)') .c.'Initial particle configuration generated'
            write (stdout, '(a,i0)') .c.'Total number of particles: ', ntotal

        end if

        close (1)
        close (2)
        close (3)

1001    format(1x, i5, 6(2x, e14.8))
1002    format(1x, i5, 7(2x, e14.8))
1003    format(1x, i5, 2x, i2, 2x, e14.8)

    end subroutine input

    !> 一维振荡管的初始数据
    subroutine shock_tube(x, vx, mass, rho, p, u, itype, hsml, ntotal)
        real(rk), intent(out) :: x(:, :)    !! 粒子位置
        real(rk), intent(out) :: vx(:, :)   !! 粒子速度
        real(rk), intent(out) :: mass(:)    !! 粒子质量
        real(rk), intent(out) :: rho(:)     !! 粒子密度
        real(rk), intent(out) :: p(:)       !! 粒子压力
        real(rk), intent(out) :: u(:)       !! 粒子内部能量
        integer, intent(out) :: itype(:)    !! 粒子类型
        real(rk), intent(out) :: hsml(:)    !! 粒子光滑长度
        integer, intent(out) :: ntotal      !! 粒子总数

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

    !> 二维剪切腔的初始数据
    subroutine shear_cavity(x, vx, mass, rho, p, u, itype, hsml, ntotal)
        real(rk), intent(out) :: x(:, :)    !! 粒子位置
        real(rk), intent(out) :: vx(:, :)   !! 粒子速度
        real(rk), intent(out) :: mass(:)    !! 粒子质量
        real(rk), intent(out) :: rho(:)     !! 粒子密度
        real(rk), intent(out) :: p(:)       !! 粒子压力
        real(rk), intent(out) :: u(:)       !! 粒子内部能量
        integer, intent(out) :: itype(:)    !! 粒子类型
        real(rk), intent(out) :: hsml(:)    !! 粒子光滑长度
        integer, intent(out) :: ntotal      !! 粒子总数

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

    !> 装载虚粒子或通过问题的几何形状产生虚粒子信息的子程序
    subroutine virt_part(itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype, keep)
        integer, intent(in) :: itimestep    !! 当前时间步长
        integer, intent(in) :: ntotal       !! 总粒子数
        integer, intent(out) :: nvirt       !! 虚拟粒子数
        real(rk), intent(inout) :: hsml(:)  !! 光滑长度
        real(rk), intent(inout) :: mass(:)  !! 粒子质量
        real(rk), intent(inout) :: x(:, :)  !! 粒子坐标
        real(rk), intent(inout) :: vx(:, :) !! 粒子速度
        real(rk), intent(inout) :: rho(:)   !! 密度
        real(rk), intent(inout) :: u(:)     !! 内部能量
        real(rk), intent(inout) :: p(:)     !! 粒子压力
        integer, intent(inout) :: itype(:)  !! 粒子类型
        logical, intent(in) :: keep         !! 是否读取存储的虚粒子信息

        integer :: i, j, d, im, mp
        real(rk) :: xl, dx, v_inf

        if (keep) then
            nvirt = saved_virt_part%nvirt
            x(:, ntotal + 1:ntotal + nvirt) = saved_virt_part%x
            vx(:, ntotal + 1:ntotal + nvirt) = saved_virt_part%vx
            mass(ntotal + 1:ntotal + nvirt) = saved_virt_part%mass
            rho(ntotal + 1:ntotal + nvirt) = saved_virt_part%rho
            p(ntotal + 1:ntotal + nvirt) = saved_virt_part%p
            hsml(ntotal + 1:ntotal + nvirt) = saved_virt_part%hsml
            u(ntotal + 1:ntotal + nvirt) = saved_virt_part%u
            itype(ntotal + 1:ntotal + nvirt) = saved_virt_part%itype
            return
        end if

        if (vp_input) then

            open (1, file=in_path//'/xv_vp.dat')
            open (2, file=in_path//'/state_vp.dat')
            open (3, file=in_path//'/other_vp.dat')
            read (1, *) nvirt
            do j = 1, nvirt
                i = ntotal + j
                read (1, *) im, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
                read (2, *) im, mass(i), rho(i), p(i), u(i)
                read (3, *) im, itype(i), hsml(i)
            end do
            close (1)
            close (2)
            close (3)

        else

            nvirt = 0
            mp = 40
            xl = 1.0e-3_rk
            dx = xl/mp
            v_inf = 1.e-3_rk

            !     monaghan type virtual particle on the upper side

            do i = 1, 2*mp + 1
                nvirt = nvirt + 1
                x(1, ntotal + nvirt) = (i - 1)*dx/2
                x(2, ntotal + nvirt) = xl
                vx(1, ntotal + nvirt) = v_inf
                vx(2, ntotal + nvirt) = 0.0_rk
            end do

            !     monaghan type virtual particle on the lower side

            do i = 1, 2*mp + 1
                nvirt = nvirt + 1
                x(1, ntotal + nvirt) = (i - 1)*dx/2
                x(2, ntotal + nvirt) = 0.0_rk
                vx(1, ntotal + nvirt) = 0.0_rk
                vx(2, ntotal + nvirt) = 0.0_rk
            end do

            !     monaghan type virtual particle on the left side

            do i = 1, 2*mp - 1
                nvirt = nvirt + 1
                x(1, ntotal + nvirt) = 0.0_rk
                x(2, ntotal + nvirt) = i*dx/2
                vx(1, ntotal + nvirt) = 0.0_rk
                vx(2, ntotal + nvirt) = 0.0_rk
            end do

            !     monaghan type virtual particle on the right side

            do i = 1, 2*mp - 1
                nvirt = nvirt + 1
                x(1, ntotal + nvirt) = xl
                x(2, ntotal + nvirt) = i*dx/2
                vx(1, ntotal + nvirt) = 0.0_rk
                vx(2, ntotal + nvirt) = 0.0_rk
            end do

            do i = 1, nvirt
                rho(ntotal + i) = 1000.0_rk
                mass(ntotal + i) = rho(ntotal + i)*dx*dx
                p(ntotal + i) = 0.0_rk
                u(ntotal + i) = 357.1_rk
                itype(ntotal + i) = -2
                hsml(ntotal + i) = dx
            end do

        end if

        if (.not. keep) then
            saved_virt_part%nvirt = nvirt
            saved_virt_part%x = x(:, ntotal + 1:ntotal + nvirt)
            saved_virt_part%vx = vx(:, ntotal + 1:ntotal + nvirt)
            saved_virt_part%mass = mass(ntotal + 1:ntotal + nvirt)
            saved_virt_part%rho = rho(ntotal + 1:ntotal + nvirt)
            saved_virt_part%hsml = hsml(ntotal + 1:ntotal + nvirt)
            saved_virt_part%p = p(ntotal + 1:ntotal + nvirt)
            saved_virt_part%u = u(ntotal + 1:ntotal + nvirt)
            saved_virt_part%itype = itype(ntotal + 1:ntotal + nvirt)
        end if

        if (mod(itimestep, save_step) == 0) then
            open (1, file=in_path//'/xv_vp.dat')
            open (2, file=in_path//'/state_vp.dat')
            open (3, file=in_path//'/other_vp.dat')
            write (1, *) nvirt
            do i = ntotal + 1, ntotal + nvirt
                write (1, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
                write (2, 1002) i, mass(i), rho(i), p(i), u(i)
                write (3, 1003) i, itype(i), hsml(i)
            end do
            close (1); close (2); close (3)
        end if

        if (mod(itimestep, print_step) == 0) then
            if (int_stat) then
                write (stdout, '(a)') .c.'Statistics: virtual boundary particles:'
                write (stdout, '(a,i0)') .c.'Number of virtual particles: ', nvirt
            end if
        end if

1001    format(1x, i6, 6(2x, e14.8))
1002    format(1x, i6, 7(2x, e14.8))
1003    format(1x, i6, 2x, i4, 2x, e14.8)

    end subroutine virt_part

end module input_m
