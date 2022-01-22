!> 用于计算具有坐标 (x) 的粒子所在的排序网格的单元格的坐标 (xgcell) 的子例程。
!>   subroutine to calculate the coordinates (xgcell) of the cell of
!>   the sorting  grid, in which the particle with coordinates (x) lies.
subroutine grid_geom(i, x, ngridx, maxgridx, mingridx, dgeomx, xgcell)

    use sph_kind, only: rk
    use parameter
    implicit none

    !> 粒子序号
    !> index of particle
    integer, intent(in) :: i
    !> 粒子的坐标
    !> coordinates of particle
    real(rk), intent(in) :: x(dim)
    !> 分割网格的个数
    !> number of sorting grid cells in x, y, z-direction
    integer, intent(in) :: ngridx(dim)
    !> 分割网格的最大坐标
    !> maximum x-, y- and z-coordinate of grid range
    real(rk), intent(in) :: maxgridx(dim)
    !> 分割网格的最小坐标
    !> minimum x-, y- and z-coordinate of grid range
    real(rk), intent(in) :: mingridx(dim)
    !> 分割网格的坐标增量
    !> x-, y- and z-expansion of grid range
    real(rk), intent(in) :: dgeomx(dim)
    !> 分割网格的坐标
    !> x-, y- and z-coordinte of sorting grid cell
    real(rk), intent(out) :: xgcell(dim)

    integer :: d

    do d = 1, 3
        xgcell(d) = 1
    end do

    do d = 1, dim
        if ((x(d) > maxgridx(d)) .or. (x(d) < mingridx(d))) then
            print *, ' >>> error <<< : particle out of range'
            print *, '    particle position: x(', i, d, ') = ', x(d)
            print *, '    range: [xmin,xmax](', d, ') =                               [', mingridx(d), ',', maxgridx(d), ']'
            stop
        else
            xgcell(d) = int(real(ngridx(d))/dgeomx(d)*(x(d) - mingridx(d)) + 1._rk)
        end if
    end do

end subroutine grid_geom
