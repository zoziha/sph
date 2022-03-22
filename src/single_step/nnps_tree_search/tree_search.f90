module tree_search_m

    use config_m, only: rk, sp, stdout
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
    subroutine tree_search(itimestep, ntotal, hsml, x, niac, pair_i, &
                           pair_j, w, dwdx, countiac)
        integer, intent(in) :: itimestep, ntotal
        real(rk), intent(in) :: hsml(maxn)
        real(rk), intent(in) :: x(dim, maxn)
        integer, intent(out) :: niac
        integer, intent(out) :: pair_i(max_interaction)
        integer, intent(out) :: pair_j(max_interaction)
        real(rk), intent(out) :: w(max_interaction)
        real(rk), intent(out) :: dwdx(dim, max_interaction)
        integer, intent(out) :: countiac(maxn)

        logical info
        real(rk), save :: range_width = 0       !! 查找域宽度
        integer scale_k, i, j, k
        real(rk) dx(dim), r !! 粒子间距
        class(shape_t), allocatable :: boundary, range !! 查找域
        save boundary !! 求解域
        type(ntree_t) :: ntree !! 树
        type(queue_t) :: found
        class(*), allocatable :: data

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        end select

        countiac(1:ntotal) = 0
        niac = 0

        ! 根据各维度最大值、最小值，计算比求解域更大的四(八)叉树边界
        if (.not. allocated(boundary)) call make_boundary(minval(real(x(:, 1:ntotal), sp), dim=2), &
                                                          maxval(real(x(:, 1:ntotal), sp), dim=2), &
                                                          tree_scale_ratio, boundary, .true.)
        call make_ntree(boundary, 1, ntree)
        do i = 1, ntotal
            call ntree%insert(point_t(x(:, i), i), info)
            if (.not.info) then
                write(stdout, '(a,i0)'), 'insert error: ', i
                error stop '*<tree_search_m::tree_search>*'
            end if
        end do

        do i = 1, ntotal - 1
            call make_range(real(x(:, i), sp), real(scale_k*hsml(i), sp), range)

            call ntree%query(range, found)
            if (found%size() == 0) cycle

            do k = 1, found%size()
                call found%dequeue(data)
                select type (data)
                type is (point_t)
                    j = data%id
                end select
                if (j <= i) cycle  ! 如果域内粒子序号小于等于当前粒子序号，说明已经计录过了，则跳过

                niac = niac + 1
                pair_i(niac) = i
                pair_j(niac) = j
                countiac(i) = countiac(i) + 1
                countiac(j) = countiac(j) + 1

                call get_distance(x(1:dim, i), x(1:dim, j), dx, r)
                !> 计算核函数值和导数备用
                call kernel(r, dx, hsml(i), w(niac), dwdx(1:dim, niac))

            end do

        end do

        !     statistics for the interaction
        call set_statistics_print(itimestep, ntotal, niac, countiac)

    end subroutine tree_search

end module tree_search_m
