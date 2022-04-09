!> 流体宏观量: 
!>
!> 质量、粒子类型、位置、速度、密度、压强、内能、光滑长度、声速、熵、总能量
module macro_m

    use config_m, only: rk, stderr
    use error_stop_m, only: error_stop
    implicit none
    private :: rk
    
    real(rk), allocatable :: x(:, :), vx(:, :), mass(:), rho(:), p(:), &
                             u(:), c(:), s(:), e(:), hsml(:)
    integer, allocatable :: itype(:)  !! 粒子类型 (1: 理想气体, 2: 淡水) @todo: 采用enum类型

contains

    !> 为流、固体宏观量分配内存
    subroutine alloc_macro_memory(dim, maxn, x, vx, mass, rho, p, u, c, s, e, hsml, itype)
        integer, intent(in) :: dim                      !! 粒子维数
        integer, intent(in) :: maxn                     !! 粒子数量
        real(rk), allocatable, intent(out) :: x(:, :)   !! 粒子位置
        real(rk), allocatable, intent(out) :: vx(:, :)  !! 粒子速度
        real(rk), allocatable, intent(out) :: mass(:)   !! 粒子质量
        real(rk), allocatable, intent(out) :: rho(:)    !! 粒子密度
        real(rk), allocatable, intent(out) :: p(:)      !! 粒子压强
        real(rk), allocatable, intent(out) :: u(:)      !! 粒子内能
        real(rk), allocatable, intent(out) :: c(:)      !! 粒子声速
        real(rk), allocatable, intent(out) :: s(:)      !! 粒子熵
        real(rk), allocatable, intent(out) :: e(:)      !! 粒子总能
        real(rk), allocatable, intent(out) :: hsml(:)   !! 粒子光滑长度
        integer, allocatable, intent(out) :: itype(:)   !! 粒子类型

        integer :: alloc_err
        allocate(x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), &
                 u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), itype(maxn), stat=alloc_err)
        if (alloc_err /= 0) then
            call error_stop('Error in allocating memory for macro arrays', &
                            'macro_m::alloc_macro_memory')
        end if
        
    end subroutine alloc_macro_memory

end module macro_m
