!>   subroutine to calculate the coordinates (xgcell) of the cell of
!>   the sorting  grid, in which the particle with coordinates (x) lies.
!>
!>     x        : coordinates of particle                            [in]
!>     ngridx   : number of sorting grid cells in x, y, z-direction  [in]
!>     maxgridx : maximum x-, y- and z-coordinate of grid range      [in]
!>     mingridx : minimum x-, y- and z-coordinate of grid range      [in]
!>     dgeomx   : x-, y- and z-expansion of grid range               [in]
!>     xgcell   : x-, y- and z-coordinte of sorting grid cell       [out]

subroutine grid_geom(i, x, ngridx, maxgridx, mingridx, dgeomx, xgcell)

    use sph_kind, only: rk
    implicit none
    include 'param.inc'

    integer  :: i, ngridx(dim), xgcell(3)
    real(rk) :: x(dim), maxgridx(dim), mingridx(dim), dgeomx(dim)
    integer  :: d

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
