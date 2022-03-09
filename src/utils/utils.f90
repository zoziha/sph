!> 通用工具
!> General Tools
module utils

    use sph_kinds, only: rk
    use parameter
    implicit none
    private

    public :: get_distance

contains

    !> 获取两点之间的距离
    subroutine get_distance(x, y, d, r)
        real(rk), intent(in), dimension(dim) :: x, y
        real(rk), intent(out), dimension(dim) :: d
        real(rk), intent(out) :: r !! 欧式距离

        d = x - y
        r = norm2(d)

    end subroutine get_distance

end module utils
