! 前处理，使用 Lua 脚本初始化粒子。
module lua_call_m

    use easy_lua_m, only: get_value
    use, intrinsic :: iso_c_binding, only: c_ptr
    use lua
    use config_m, only: rk, in_path
    use input_m, only: saved_virt_part
    use parameter, only: dim
    implicit none
    private

    public :: lua_input, lua_virt_part

contains

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

    ! 读取虚粒子数据
    subroutine lua_virt_part(x, vx, mass, rho, p, u, itype, hsml, ntotal, nvirt, keep)
        use config_m, only: lua_script
        integer, intent(in) :: ntotal
        !> 虚拟粒子数
        !> Number of virtual particles
        integer, intent(out) :: nvirt
        !> 光滑长度
        !> Smoothing length
        real(rk), intent(inout) :: hsml(:)
        !> 粒子质量
        !> Particle masses
        real(rk), intent(inout) :: mass(:)
        !> 粒子坐标
        !> Particle coordinates
        real(rk), intent(inout) :: x(:, :)
        !> 粒子速度
        !> Particle velocities
        real(rk), intent(inout) :: vx(:, :)
        !> 密度
        !> Density
        real(rk), intent(inout) :: rho(:)
        !> 内部能量
        !> Internal energy
        real(rk), intent(inout) :: u(:)
        !> 粒子压力
        !> Particle pressure
        real(rk), intent(inout) :: p(:)
        !> 粒子类型
        !> Particle type
        integer, intent(inout) :: itype(:)
        logical, intent(in) :: keep
        type(c_ptr) :: l
        integer :: rc

        if (keep) then
            nvirt = saved_virt_part%nvirt
            x(:, ntotal + 1:ntotal + nvirt) = saved_virt_part%x
            vx(:, ntotal + 1:ntotal + nvirt) = saved_virt_part%vx
            mass(ntotal + 1:ntotal + nvirt) = saved_virt_part%mass
            rho(ntotal + 1:ntotal + nvirt) = saved_virt_part%rho
            p(ntotal + 1:ntotal + nvirt) = saved_virt_part%p
            hsml(ntotal + 1:ntotal + nvirt) = saved_virt_part%hsml
            u(ntotal + 1:ntotal + nvirt) = saved_virt_part%u
            itype(ntotal + 1:ntotal + nvirt) = saved_virt_part%itype
            return
        else
            l = lual_newstate()
            call lual_openlibs(l)
            rc = lual_loadfile(l, in_path//'/'//lua_script)
            rc = lua_pcall(l, 0, 0, 0)
            rc = lua_getglobal(l, "virt_part")
            rc = lua_pcall(l, 0, 9, 0)
            call get_value(l, saved_virt_part%x, index=-9) ! 倒序读取 Lua 堆栈
            call get_value(l, saved_virt_part%vx, index=-8)
            call get_value(l, saved_virt_part%mass, index=-7)
            call get_value(l, saved_virt_part%rho, index=-6)
            call get_value(l, saved_virt_part%p, index=-5)
            call get_value(l, saved_virt_part%u, index=-4)
            call get_value(l, saved_virt_part%itype, index=-3)
            call get_value(l, saved_virt_part%hsml, index=-2)
            call get_value(l, saved_virt_part%nvirt, index=-1)
            call lua_close(l)
        end if

    end subroutine lua_virt_part

end module lua_call_m
