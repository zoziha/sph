!> 保存粒子信息到硬盘文件中的子程序
!> subroutine for saving particle information to external disk file
subroutine output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)

    use sph_kind, only: rk
    use parameter
    implicit none

    !> 粒子的位置
    !> coordinates of particles
    real(rk), intent(in) :: x(dim, maxn)
    !> 粒子的速度
    !> velocities of particles
    real(rk), intent(in) :: vx(dim, maxn)
    !> 粒子的质量
    !> mass of particles
    real(rk), intent(in) :: mass(maxn)
    !> 粒子的密度
    !> dnesities of particles
    real(rk), intent(in) :: rho(maxn)
    !> 粒子的压力
    !> pressure  of particles
    real(rk), intent(in) :: p(maxn)
    !> 粒子的内部能量
    !> internal energy of particles
    real(rk), intent(in) :: u(maxn)
    !> 粒子的声速
    !> sound velocity of particles
    real(rk), intent(in) :: c(maxn)
    !> 粒子的类型 (1: ideal gas; 2: water; 3: TNT)
    !> types of particles
    integer, intent(in) :: itype(maxn)
    !> 粒子的光滑长度
    !> smoothing lengths of particles
    real(rk), intent(in) :: hsml(maxn)
    !> 粒子的总数
    !> total particle number
    integer, intent(in) :: ntotal

    integer :: i, d, npart

    open (1, file='./data/f_xv.dat')
    open (2, file='./data/f_state.dat')
    open (3, file='./data/f_other.dat')

    write (1, *) ntotal
    do i = 1, ntotal
        write (1, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
        write (2, 1002) i, mass(i), rho(i), p(i), u(i)
        write (3, 1003) i, itype(i), hsml(i)
    end do

    close (1)
    close (2)
    close (3)

1001 format(1x, i6, 6(2x, e14.8))
1002 format(1x, i6, 7(2x, e14.8))
1003 format(1x, i6, 2x, i4, 2x, e14.8)

end subroutine output
