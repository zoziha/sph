!> 装载虚粒子或通过问题的几何形状产生虚粒子信息的子程序。
!> @todo: 不重复输出虚粒子信息
!>   subroutine to determine the information of virtual particles
!>   here only the monaghan type virtual particles for the 2d shear
!>   cavity driven problem are generated.
subroutine virt_part(itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype)

    use config_m, only: rk, stdout
    use info_m, only: operator(.c.)
    use parameter
    implicit none

    !> 当前时间步
    !> Current time step
    integer, intent(in) :: itimestep
    !> 总粒子数
    !> Total number of particles
    integer, intent(in) :: ntotal
    !> 虚拟粒子数
    !> Number of virtual particles
    integer, intent(out) :: nvirt
    !> 光滑长度
    !> Smoothing length
    real(rk), intent(inout) :: hsml(maxn)
    !> 粒子质量
    !> Particle masses
    real(rk), intent(inout) :: mass(maxn)
    !> 粒子坐标
    !> Particle coordinates
    real(rk), intent(inout) :: x(dim, maxn)
    !> 粒子速度
    !> Particle velocities
    real(rk), intent(inout) :: vx(dim, maxn)
    !> 密度
    !> Density
    real(rk), intent(inout) :: rho(maxn)
    !> 内部能量
    !> Internal energy
    real(rk), intent(inout) :: u(maxn)
    !> 粒子压力
    !> Particle pressure
    real(rk), intent(inout) :: p(maxn)
    !> 粒子类型
    !> Particle type
    integer, intent(inout) :: itype(maxn)

    integer :: i, j, d, im, mp
    real(rk) :: xl, dx, v_inf

    if (vp_input) then

        open (1, file='./data/xv_vp.dat')
        open (2, file='./data/state_vp.dat')
        open (3, file='./data/other_vp.dat')
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

    if (mod(itimestep, save_step) == 0) then
        open (1, file='./data/xv_vp.dat')
        open (2, file='./data/state_vp.dat')
        open (3, file='./data/other_vp.dat')
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
            write(stdout, '(a)') .c.'Statistics: virtual boundary particles:'
            write(stdout, '(a,i0)')  .c.'Number of virtual particles: ', nvirt
        end if
    end if

1001 format(1x, i6, 6(2x, e14.8))
1002 format(1x, i6, 7(2x, e14.8))
1003 format(1x, i6, 2x, i4, 2x, e14.8)

end subroutine virt_part
