! 流体宏观量: 质量、粒子类型、位置、速度、密度、压强、内能、光滑长度、声速、熵、总能量
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
        implicit none
        integer, intent(in) :: dim, maxn
        real(rk), allocatable, intent(out) :: x(:, :), vx(:, :), mass(:), rho(:), p(:), &
                                             u(:), c(:), s(:), e(:), hsml(:)
        integer, allocatable, intent(out) :: itype(:)
        integer :: alloc_err
        allocate(x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), &
                 u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), itype(maxn), stat=alloc_err)
        if (alloc_err /= 0) then
            call error_stop('Error in allocating memory for macro arrays', &
                            'macro_m::alloc_macro_memory')
        end if
        
    end subroutine alloc_macro_memory

end module macro_m
