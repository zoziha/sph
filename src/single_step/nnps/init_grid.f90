!>   subroutine to established a pair linked list by sorting grid cell.
!>   it is suitable for a homogeneous particle distribution with the
!>   same smoothing length in an instant. a fixed number of particles
!>   lie in each cell.
!>
!>     ntotal   : number of particles                                [in]
!>     hsml     : smoothing length                                   [in]
!>     grid     : array of grid cells                               [out]
!>     ngridx   : number of sorting grid cells in x, y, z-direction [out]
!>     ghsmlx   : smoothing length measured in cells of the grid    [out]
!>     maxgridx : maximum x-, y- and z-coordinate of grid range     [out]
!>     mingridx : minimum x-, y- and z-coordinate of grid range     [out]
!>     dgeomx   : x-, y- and z-expansion of grid range              [out]

subroutine init_grid(ntotal, hsml, grid, ngridx, ghsmlx, maxgridx, mingridx, dgeomx)

    use sph_kind, only: rk
    use parameter
    implicit none

    !     parameter used for sorting grid cells in the link list algorithm
    !     maxngx  : maximum number of sorting grid cells in x-direction
    !     maxngy  : maximum number of sorting grid cells in y-direction
    !     maxngz  : maximum number of sorting grid cells in z-direction
    !     determining maximum number of sorting grid cells:
    !     (for an homogeneous particle distribution:)
    !     1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
    !     2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
    !     3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)
    integer :: maxngx, maxngy, maxngz
    parameter(maxngx=100, maxngy=100, maxngz=1)
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    integer :: grid(maxngx, maxngy, maxngz), ngridx(dim), ghsmlx(dim)
    real(rk) :: hsml, maxgridx(dim), mingridx(dim), dgeomx(dim)
    integer :: i, j, k, d, maxng(dim), ngrid(3)
    !> averaged number of particles per grid cell
    real(rk), parameter :: nppg = 3.0_rk

    !     initialize parameters: maximum number of grid cells

    maxng(1) = maxngx
    if (dim >= 2) then
        maxng(2) = maxngy
        if (dim == 3) then
            maxng(3) = maxngz
        end if
    end if

    do d = 1, 3
        ngrid(d) = 1
    end do

    !     range of sorting grid

    maxgridx(1) = x_maxgeom
    mingridx(1) = x_mingeom
    if (dim >= 2) then
        maxgridx(2) = y_maxgeom
        mingridx(2) = y_mingeom
        if (dim == 3) then
            maxgridx(3) = z_maxgeom
            mingridx(3) = z_mingeom
        end if
    end if

    do d = 1, dim
        dgeomx(d) = maxgridx(d) - mingridx(d)
    end do

    !     number of grid cells in x-, y- and z-direction:

    if (dim == 1) then
        ngridx(1) = min(int(ntotal/nppg) + 1, maxng(1))
    else if (dim == 2) then
        ngridx(1) = min(int(sqrt(ntotal*dgeomx(1)/(dgeomx(2)*nppg))) + 1, maxng(1))
        ngridx(2) = min(int(ngridx(1)*dgeomx(2)/dgeomx(1)) + 1, maxng(2))
    else if (dim == 3) then
        ngridx(1) = min(int((ntotal*dgeomx(1)*dgeomx(1)/(dgeomx(2)*dgeomx(3)*nppg))**(1._rk/3._rk)) + 1, maxng(1))
        ngridx(2) = min(int(ngridx(1)*dgeomx(2)/dgeomx(1)) + 1, maxng(2))
        ngridx(3) = min(int(ngridx(1)*dgeomx(3)/dgeomx(1)) + 1, maxng(3))
    end if

    !     smoothing length measured in grid cells:

    do d = 1, dim
        ghsmlx(d) = int(real(ngridx(d))*hsml/dgeomx(d)) + 1
    end do

    do d = 1, dim
        ngrid(d) = ngridx(d)
    end do

    !     initialize grid

    do i = 1, ngrid(1)
        do j = 1, ngrid(2)
            do k = 1, ngrid(3)
                grid(i, j, k) = 0
            end do
        end do
    end do

end subroutine init_grid
