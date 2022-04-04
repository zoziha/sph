module link_list_m

    use config_m, only: rk, skf, max_interaction
    use parameter, only: dim
    use output_m, only: set_statistics_print
    use kernel_m, only: kernel
    implicit none
    private

    public :: link_list

contains

    !>   subroutine to calculate the smoothing funciton for each particle and
    !>   the interaction parameters used by the sph algorithm. interaction
    !>   pairs are determined by using a sorting grid linked list
    !>
    !>     itimestep : current time step                                 [in]
    !>     ntotal    : number of particles                               [in]
    !>     hsml      : smoothing length, same for all particles          [in]
    !>     x         : coordinates of all particles                      [in]
    !>     niac      : number of interaction pairs                      [out]
    !>     pair_i    : list of first partner of interaction pair        [out]
    !>     pair_j    : list of second partner of interaction pair       [out]
    !>     w         : kernel for all interaction pairs                 [out]
    !>     dwdx      : derivative of kernel with respect to x, y and z  [out]
    !>     countiac  : number of neighboring particles                  [out]
    subroutine link_list(itimestep, ntotal, hsml, x, niac, pair_i, pair_j, w, dwdx, countiac)
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
        integer, intent(in) :: itimestep
        !> 在模拟中所使用的粒子总数
        !> number of particles in simulation
        integer, intent(in) :: ntotal
        !> 相互作用对的数目
        integer, intent(out) :: niac
        integer :: pair_i(:), pair_j(:), countiac(:)
        real(rk) :: hsml, x(:, :), w(:), dwdx(:, :)
        integer :: i, j, d, scale_k
        integer :: grid(maxngx, maxngy, maxngz), xgcell(3, ntotal), gcell(3), xcell, ycell, zcell, celldata(ntotal), minxcell(3), &
                   maxxcell(3), dnxgcell(dim), dpxgcell(dim), ngridx(dim), ghsmlx(dim)
        real(rk) :: dr, r, dx(dim), mingridx(dim), maxgridx(dim), tdwdx(dim), dgeomx(dim)

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        countiac(1:ntotal) = 0

        !     initialize grid:

        call init_grid(ntotal, hsml, grid, ngridx, ghsmlx, maxgridx, mingridx, dgeomx)

        !     position particles on grid and create linked list:

        do i = 1, ntotal
            call grid_geom(i, x(:, i), ngridx, maxgridx, mingridx, dgeomx, gcell)
            do d = 1, dim
                xgcell(d, i) = gcell(d)
            end do
            celldata(i) = grid(gcell(1), gcell(2), gcell(3))
            grid(gcell(1), gcell(2), gcell(3)) = i
        end do

        !     determine interaction parameters:

        niac = 0
        do i = 1, ntotal - 1

            !     determine range of grid to go through:

            do d = 1, 3
                minxcell(d) = 1
                maxxcell(d) = 1
            end do
            do d = 1, dim
                dnxgcell(d) = xgcell(d, i) - ghsmlx(d)
                dpxgcell(d) = xgcell(d, i) + ghsmlx(d)
                minxcell(d) = max(dnxgcell(d), 1)
                maxxcell(d) = min(dpxgcell(d), ngridx(d))
            end do

            !     search grid:

            do zcell = minxcell(3), maxxcell(3)
                do ycell = minxcell(2), maxxcell(2)
                    do xcell = minxcell(1), maxxcell(1)
                        j = grid(xcell, ycell, zcell)
1                       if (j > i) then
                            dx(1) = x(1, i) - x(1, j)
                            dr = dx(1)*dx(1)
                            do d = 2, dim
                                dx(d) = x(d, i) - x(d, j)
                                dr = dr + dx(d)*dx(d)
                            end do
                            if (sqrt(dr) < scale_k*hsml) then
                                if (niac < max_interaction) then

                                    !     neighboring pair list, and totalinteraction number and
                                    !     the interaction number for each particle

                                    niac = niac + 1
                                    pair_i(niac) = i
                                    pair_j(niac) = j
                                    r = sqrt(dr)
                                    countiac(i) = countiac(i) + 1
                                    countiac(j) = countiac(j) + 1

                                    !--- kernel and derivations of kernel

                                    call kernel(r, dx, hsml, w(niac), tdwdx)
                                    do d = 1, dim
                                        dwdx(d, niac) = tdwdx(d)
                                    end do
                                else
                                    error stop ' >>> error <<< : too many interactions'
                                end if
                            end if
                            j = celldata(j)
                            goto 1
                        end if
                    end do
                end do
            end do
        end do

        !     statistics for the interaction
        call set_statistics_print(itimestep, ntotal, niac, countiac)

    end subroutine link_list

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

        use config_m, only: rk
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
        real(rk) :: hsml, maxgridx(:), mingridx(:), dgeomx(:)
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
!> 用于计算具有坐标 (x) 的粒子所在的排序网格的单元格的坐标 (xgcell) 的子例程。
!>   subroutine to calculate the coordinates (xgcell) of the cell of
!>   the sorting  grid, in which the particle with coordinates (x) lies.
    subroutine grid_geom(i, x, ngridx, maxgridx, mingridx, dgeomx, xgcell)

        use config_m, only: rk
        use parameter
        implicit none

        !> 粒子序号
        !> index of particle
        integer, intent(in) :: i
        !> 粒子的坐标
        !> coordinates of particle
        real(rk), intent(in) :: x(:)
        !> 分割网格的个数
        !> number of sorting grid cells in x, y, z-direction
        integer, intent(in) :: ngridx(:)
        !> 分割网格的最大坐标
        !> maximum x-, y- and z-coordinate of grid range
        real(rk), intent(in) :: maxgridx(:)
        !> 分割网格的最小坐标
        !> minimum x-, y- and z-coordinate of grid range
        real(rk), intent(in) :: mingridx(:)
        !> 分割网格的坐标增量
        !> x-, y- and z-expansion of grid range
        real(rk), intent(in) :: dgeomx(:)
        !> 分割网格的坐标
        !> x-, y- and z-coordinte of sorting grid cell
        integer, intent(out) :: xgcell(3)

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
                xgcell(d) = int(real(ngridx(d))/dgeomx(d)*(x(d) - mingridx(d)) + 1.0_rk)
            end if
        end do

    end subroutine grid_geom

end module link_list_m
