module tree_search_m

    use config_m, only: rk, stdout, tinsert, tsearch, skf
    use queue_m, only: queue_t
    use parameter
    use output_m, only: set_statistics_print
    use utils, only: get_distance
    use ntree_factory_m, only: ntree_t, shape_t, point_t, &
                               make_ntree, make_boundary, make_range
    implicit none
    private

    public :: tree_search

contains

    !> 树型搜索法，适用于变光滑长度 @todo: 3维
    !@todo: 使用记忆<n>叉树，提高搜索速度
    subroutine tree_search(itimestep, ntotal, hsml, x, niac, pair_i, &
                           pair_j, w, dwdx, countiac)
        integer, intent(in) :: itimestep, ntotal
        real(rk), intent(in) :: hsml(:)
        real(rk), intent(in) :: x(:, :)
        integer, intent(out) :: niac
        integer, intent(out) :: pair_i(:)
        integer, intent(out) :: pair_j(:)
        real(rk), intent(out) :: w(:)
        real(rk), intent(out) :: dwdx(:, :)
        integer, intent(out) :: countiac(:)

        logical info
        integer scale_k, i, k
        real(rk) dx(dim), r !! 粒子间距
        class(shape_t), allocatable :: boundary, range !! 查找域
        save boundary !! 求解域
        type(ntree_t) :: ntree !! 树
        type(queue_t) :: found
        class(*), allocatable :: j
        real :: t1, t2

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        countiac(1:ntotal) = 0
        niac = 0

        ! 根据各维度最大值、最小值，计算比求解域更大的四 (八) 叉树边界
        if (.not. allocated(boundary)) call make_boundary(minval(x(:, 1:ntotal), dim=2), &
                                                          maxval(x(:, 1:ntotal), dim=2), &
                                                          tree_scale_ratio, boundary, .true.)
        call make_ntree(boundary, 4, ntree)  !! 搜索花费时间较长，可以考虑减小树的深度，设置为2~4
        call cpu_time(t1)
        do i = 1, ntotal
            call ntree%insert(point_t(x=x(1:dim, i), id=i), info)
            if (.not. info) then
                write (stdout, '(a,i0)') 'insert error: ', i
                error stop '*<tree_search_m::tree_search>*'
            end if
        end do
        call cpu_time(t2)
        tinsert = t2 - t1 + tinsert

        do i = 1, ntotal - 1
            call make_range(x(:, i), scale_k*hsml(i), range)

            call ntree%query(range, found)
            if (found%size() == 0) cycle

            do k = 1, found%size()
                call found%dequeue(j)
                select type (j)
                type is (integer)
                    if (j <= i) cycle  ! 如果域内粒子序号小于等于当前粒子序号，说明已经计录过了，则跳过
                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1

                    call get_distance(x(1:dim, i), x(1:dim, j), dx, r)
                    !> 计算核函数值和导数备用
                    call kernel(r, dx, hsml(i), w(niac), dwdx(1:dim, niac))
                class default
                    write (stdout, '(a,i0)') 'error item type: ', i
                    error stop '*<tree_search_m::tree_search>*'
                end select

            end do

        end do
        call cpu_time(t1)
        tsearch = t1 - t2 + tsearch

        !     statistics for the interaction
        call set_statistics_print(itimestep, ntotal, niac, countiac)

    end subroutine tree_search

end module tree_search_m
