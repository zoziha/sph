module test_tree_search_m

    use testdrive, only: new_unittest, unittest_type, error_type, check
    use tree_search_m, only: tree_search
    use quad_types, only: point_t, rectangle_t, circle_t
    use quad, only: quad_tree_t
    use sph_kinds, only: rk
    implicit none
    private

    public :: collect_tree_search_tests

    type(quad_tree_t) :: qt

contains

    subroutine collect_tree_search_tests(test_suite)
        type(unittest_type), allocatable, intent(out) :: test_suite(:)

        test_suite = [ &
                     new_unittest("test_set_quad_tree", test_set_quad_tree), &
                     new_unittest("test_get_query_tree", test_get_query_tree) &
                     ]

    end subroutine collect_tree_search_tests

    subroutine test_set_quad_tree(error)
        type(error_type), allocatable, intent(out) :: error

        call qt%constructor(rectangle_t(0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk), 1)
        call check(error, qt%insert(point_t(0.5_rk, 0.5_rk, 1)), "insert 1")
        if (allocated(error)) return

        call check(error, qt%insert(point_t(-0.3_rk, 0.7_rk, 2)), "insert 2")
        if (allocated(error)) return

        call check(error, qt%insert(point_t(0.7_rk, -0.3_rk, 3)), "insert 3")
        if (allocated(error)) return

        call check(error, qt%insert(point_t(-0.1_rk, 0.999_rk, 4)), "insert 4")
        if (allocated(error)) return

        call check(error, qt%insert(point_t(0.0_rk, -0.0_rk, 5)), "insert 5")
        if (allocated(error)) return
        
        call check(error, qt%insert(point_t(0.0_rk, 0.5_rk, 6)), "insert 6")
        
        !> 域外点
        call check(error, .not.qt%insert(point_t(1.0_rk, 2.0_rk, 6)), "insert 6")

    end subroutine test_set_quad_tree

    subroutine test_get_query_tree(error)
        type(error_type), allocatable, intent(out) :: error
        type(circle_t) range
        type(point_t), allocatable :: found(:)
        integer is, index(5)
        
        index = [1, 5, 6, 2, 3]
        range = circle_t(0.0_rk, 0.0_rk, 1.0_rk)
        call qt%query(range, found)
        call check(error, size(found), 5, "query size")
        if (allocated(error)) return
        
        do is = 1, size(found)
            call check(error, found(is)%index, index(is), "check index")
            if (allocated(error)) return
        end do
        deallocate(found)
        
        range = circle_t(0.1_rk, -0.2_rk, 0.6_rk)
        call qt%query(range, found)
        call check(error, size(found), 1, "query size")
        if (allocated(error)) return
        
        call check(error, found(1)%index, 5, "check index")

    end subroutine test_get_query_tree

end module test_tree_search_m
