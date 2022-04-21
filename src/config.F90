!> SPH 精度
module config_m

    use, intrinsic :: iso_fortran_env, only: &
        sp => real32, dp => real64, stdout => output_unit, &
        stdin => input_unit, stderr => error_unit
    implicit none
    public

    character(*), parameter :: version = '0.2.0'    !! sph solver 版本
    integer, parameter :: rk = dp           !! sph预设浮点精度 @todo: not compatible with single precision at this time

    real, save :: tinsert = 0.0             !! 粒子记录时间
    real, save :: tsearch = 0.0             !! 粒子搜索时间

    character(:), allocatable :: in_path    !! 输入路径
    character(:), allocatable :: out_path   !! 输出路径

    character(:), allocatable :: nick       !! 工程名
    real(rk) :: dt                          !! 时间步长
    real(rk) :: CFL                         !! CFL条件数值, 默认值 0.3
    integer :: skf                          !! 光滑核函数的指示变量
    integer :: nnps                         !! 最近相邻粒子搜索算法指示变量

    integer :: print_step                   !! 输出到屏幕的时间步间隔
    integer :: save_step                    !! 输出到磁盘的时间步间隔
    integer :: moni_particle                !! 监控粒子, 默认值 1600

    integer :: maxn                         !! 粒子总数
    integer :: max_interaction              !! 粒子最大互动数
    logical :: virtual_part                 !! 是否使用虚拟粒子, 默认值 T
    logical :: self_gravity                 !! 是否考虑自重, 默认值 F
    logical :: visc                         !! 是否考虑粘性, 默认值 T
    integer :: eos_form                     !! 水的 EOS 形式, 默认值 2
    real(rk) :: B                           !! 弱可压 EOS 参数, 默认值 ?
    real(rk) :: rho0                        !! 参考密度, 默认值 1000.0
    real(rk) :: h_SWL                       !! 静水面高度, 默认值 1.0
    real(rk) :: c                           !! 人工声速, 默认值 nil
    logical :: visc_artificial              !! 是否考虑人工粘性, 默认值 F
    logical :: heat_artificial              !! 是否考虑人工热量, 默认值 F

    logical :: dofile                       !! 是否使用脚本, 默认值 F
    character(:), allocatable :: lua_script !! 脚本名

end module config_m
