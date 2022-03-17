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

    use config_m, only: rk
    use parameter
    use output_m, only: set_statistics_print
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
    integer, intent(in) :: itimestep
    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 相互作用对的数目
    integer, intent(out) :: niac
    integer :: pair_i(max_interaction), pair_j(max_interaction), countiac(maxn)
    real(rk) :: hsml, x(dim, maxn), w(max_interaction), dwdx(dim, max_interaction)
    integer :: i, j, d, scale_k
    integer :: grid(maxngx, maxngy, maxngz), xgcell(3, maxn), gcell(3), xcell, ycell, zcell, celldata(maxn), minxcell(3), &
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
        call grid_geom(i, x(1, i), ngridx, maxgridx, mingridx, dgeomx, gcell)
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
1                   if (j > i) then
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
