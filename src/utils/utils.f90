!> 通用工具
module utils

    use config_m, only: rk
    use parameter
    implicit none
    private

    public :: get_distance

contains

    !> 获取两点之间的距离
    pure subroutine get_distance(x, y, d, r)
        real(rk), intent(in), dimension(:) :: x, y  !! 两点坐标
        real(rk), intent(out), dimension(:) :: d    !! 坐标轴距离
        real(rk), intent(out) :: r                  !! 欧式距离

        d = x - y
        r = norm2(d)

    end subroutine get_distance

end module utils
