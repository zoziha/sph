module physical_constant_m

    use config_m, only: rk
    implicit none
    private :: rk
    
    real(rk), parameter :: g = 9.8_rk  !! 重力加速度

end module physical_constant_m