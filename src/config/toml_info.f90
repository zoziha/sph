! 读取toml配置文件
module toml_info_m

    use tomlf, only: toml_table, get_value, toml_parse
    use swift_file_m, only: is_exist
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use config_m
    use error_stop_m, only: error_stop
    implicit none
    private

    public :: parse_toml_info

contains

    !> 解析程序所需的所有配置文件
    subroutine parse_toml_info()
        call parse_access_toml()
        call parse_info_sph()
    end subroutine parse_toml_info

    ! 解析access.toml文件
    subroutine parse_access_toml()
        type(toml_table), allocatable :: access_table
        type(toml_table), pointer :: subtable
        integer access_toml_unit
        if (.not. is_exist("access.toml")) &
            call error_stop('access.toml文件不存在, 请检查!', &
                            'toml_info_m%parse_access_toml')
        open (newunit=access_toml_unit, file='access.toml', status='old')
        call toml_parse(access_table, access_toml_unit)
        close (access_toml_unit)

        call get_value(access_table, 'access', subtable)
        call get_value(subtable, 'in', in_path, './data/input')  !! 默认输入路径
        call get_value(subtable, 'out', out_path, './data/output')  !! 默认输出路径

        nullify (subtable)
    end subroutine parse_access_toml

    !> 解析sph.toml配置文件
    subroutine parse_info_sph()
        type(toml_table), allocatable :: sph_table
        type(toml_table), pointer :: subtable
        integer :: kpair
        integer sph_toml_unit
        if (.not. is_exist(in_path//"/sph.toml")) &
            call error_stop('sph.toml文件不存在, 请检查!', &
                            'toml_info_m%parse_info_sph')
        open (newunit=sph_toml_unit, file=in_path//"/sph.toml", status='old')
        call toml_parse(sph_table, sph_toml_unit)
        close (sph_toml_unit)

        call get_value(sph_table, 'name', nick, 'untitled')  !! 默认项目名
        call get_value(sph_table, 'parameter', subtable)
        call get_value(subtable, 'dt', dt) !@todo: 设置默认值或提醒
        call get_value(subtable, 'skf', skf, 1)
        call get_value(subtable, 'nnps', nnps, 1) ! 默认直接搜索
        call get_value(subtable, 'print_step', print_step, 100)
        call get_value(subtable, 'save_step', save_step, 500)
        call get_value(subtable, 'maxn', maxn, 5000)
        call get_value(subtable, 'kpair', kpair, 20)
        max_interaction = kpair*maxn

        nullify (subtable)
    end subroutine parse_info_sph

end module toml_info_m
