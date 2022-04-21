!> 通用工具
module utils

    use config_m, only: rk
    use parameter
    use physical_constant_m, only: g
    implicit none
    private

    public :: get_distance, get_speed_of_voice

contains

    !> 获取两点之间的距离
    pure subroutine get_distance(x, y, d, r)
        real(rk), intent(in), dimension(:) :: x, y  !! 两点坐标
        real(rk), intent(out), dimension(:) :: d    !! 坐标轴距离
        real(rk), intent(out) :: r                  !! 欧式距离

        d = x - y
        r = norm2(d)

    end subroutine get_distance
    
    !> 获取人工声速
    pure subroutine get_speed_of_voice(h, rho, c)
        real(rk), intent(in) :: h   !! 静水高度
        real(rk), intent(in) :: rho !! 密度
        real(rk), intent(out) :: c  !! 声速
        real(rk) :: uMax, pMax
        
        uMax = sqrt(2 * g * h)
        pMax = rho * g * h    !@todo: 可能需要修改
        
        c = 10*max(uMax, sqrt(pMax / rho))
    end subroutine get_speed_of_voice

end module utils
