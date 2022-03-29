!> 通用工具
!> General Tools
module utils

    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: get_distance

contains

    !> 获取两点之间的距离
    pure subroutine get_distance(x, y, d, r)
        real(rk), intent(in), dimension(:) :: x, y
        real(rk), intent(out), dimension(:) :: d
        real(rk), intent(out) :: r !! 欧式距离

        d = x - y
        r = norm2(d)

    end subroutine get_distance

end module utils
