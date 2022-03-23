module test_tree_search_m

    use testdrive, only: new_unittest, unittest_type, error_type, check
    use tree_search_m, only: tree_search
    use shape_m, only: point_t, rectangle_t, circle_t
    use ntree_m, only: ntree_t
    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: collect_tree_search_tests

contains

    subroutine collect_tree_search_tests(test_suite)
        type(unittest_type), allocatable, intent(out) :: test_suite(:)

        test_suite = [ &
                     new_unittest("test_tree_search_func", test_tree_search_func) &
                     ]

    end subroutine collect_tree_search_tests

    subroutine test_tree_search_func(error)
        type(error_type), allocatable, intent(out) :: error
        real(rk) :: hsml(10), x(dim, 10), w(100), dwdx(dim, 100)
        integer :: niac, ns(10), pair_i(100), pair_j(100), ns_(10), i
        ! Ready
        x = reshape([real(rk) :: 0.18459153E-04, 0.18459158E-04, &
                     0.25928249E-04, 0.41953272E-04, &
                     0.26889891E-04, 0.65649204E-04, &
                     0.27437717E-04, 0.89790037E-04, &
                     0.27706048E-04, 0.11427938E-03, &
                     0.27869843E-04, 0.13893146E-03, &
                     0.27982794E-04, 0.16367276E-03, &
                     0.28065106E-04, 0.18847140E-03, &
                     0.28128031E-04, 0.21330950E-03, &
                     0.28177749E-04, 0.23817630E-03], [2, 10])
        hsml = 0.25000000E-04_rk
        ns_ = [2, 3, 4, 4, 4, 4, 4, 4, 3, 2] ! 第四个
        call tree_search(1, 10, hsml, x, niac, pair_i, pair_j, w, dwdx, ns)
        call check(error, niac, 17, "tree_search: niac")
        if (allocated(error)) return
        do i = 1, size(ns)
            call check(error, ns(i), ns_(i), "tree_search: ns")
            if (allocated(error)) return
        end do

    end subroutine test_tree_search_func

end module test_tree_search_m
