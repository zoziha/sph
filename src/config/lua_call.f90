! 前处理，使用 Lua 脚本初始化粒子。
module lua_call_m

    use easy_lua_m, only: get_value
    use, intrinsic :: iso_c_binding, only: c_ptr
    use lua
    use config_m, only: rk, in_path
    use parameter, only: dim
    implicit none
    private

    public :: lua_input!, lua_virt_part

contains

    !@todo: lua_input, lua_nvirt?
    !> 读取初始化粒子数据
    subroutine lua_input(lua_script, x, vx, mass, rho, p, u, itype, hsml, ntotal)
        character(*), intent(in) :: lua_script
        real(rk), intent(out) :: x(:, :), vx(:, :), mass(:), rho(:), p(:), u(:), hsml(:)
        integer, intent(out) :: itype(:), ntotal
        real(rk), allocatable :: x_(:, :), vx_(:, :), mass_(:), rho_(:), p_(:), u_(:), hsml_(:)
        integer, allocatable :: itype_(:)
        integer ntotal_
        type(c_ptr) :: l
        integer :: rc

        l = lual_newstate()
        call lual_openlibs(l)
        rc = lual_loadfile(l, in_path//'/'//lua_script)
        rc = lua_pcall(l, 0, 0, 0)
        rc = lua_getglobal(l, "input")
        rc = lua_pcall(l, 0, 9, 0)
        call get_value(l, x_, index=-9) ! 倒序读取 Lua 堆栈
        call get_value(l, vx_, index=-8)
        call get_value(l, mass_, index=-7)
        call get_value(l, rho_, index=-6)
        call get_value(l, p_, index=-5)
        call get_value(l, u_, index=-4)
        call get_value(l, itype_, index=-3)
        call get_value(l, hsml_, index=-2)
        call get_value(l, ntotal_, index=-1)
        call lua_close(l)

        !@todo: 是否需要优化
        x(1:dim, 1:ntotal_) = x_
        vx(1:dim, 1:ntotal_) = vx_
        mass(1:ntotal_) = mass_
        rho(1:ntotal_) = rho_
        p(1:ntotal_) = p_
        u(1:ntotal_) = u_
        itype(1:ntotal_) = itype_
        hsml(1:ntotal_) = hsml_
        ntotal = ntotal_

    end subroutine lua_input

    ! !@todo: 处理好缝合数据
    ! subroutine lua_virt_part(lua_script, x, vx, mass, rho, p, u, itype, hsml, ntotal, nvirt)
    !     character(*), intent(in) :: lua_script
    !     real(rk), intent(inout) :: x(:, :), vx(:, :), mass(:), rho(:), p(:), u(:), hsml(:)
    !     integer, intent(inout) :: itype(:)
    !     integer, intent(in) :: ntotal
    !     integer, intent(out) :: nvirt
    !     type(c_ptr) :: l
    !     integer :: rc

    !     l = lual_newstate()
    !     call lual_openlibs(l)
    !     rc = lual_loadfile(l, lua_script)
    !     rc = lua_pcall(l, 0, 0, 0)
    !     rc = lua_getglobal(l, "virt_part")
    !     rc = lua_pcall(l, 0, 9, 0)
    !     call get_value(l, x, -1)
    !     call get_value(l, vx, -2)
    !     call get_value(l, mass, -3)
    !     call get_value(l, rho, -4)
    !     call get_value(l, p, -5)
    !     call get_value(l, u, -6)
    !     call get_value(l, itype, -7)
    !     call get_value(l, hsml, -8)
    !     call get_value(l, ntotal, -9)
    !     call get_value(l, nvirt, -10)
    !     call lual_close(l)
    ! end subroutine lua_virt_part

end module lua_call_m
