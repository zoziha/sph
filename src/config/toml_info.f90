!> 读取toml配置文件
module toml_info_m

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use arg_m, only: get_path
    use config_m
    use easy_string_m, only: to_string
    use error_stop_m, only: error_stop
    use stdlib_logger, only: stdlog => global_logger
    use swift_file_m, only: is_exist
    use tomlf, only: toml_table, get_value, toml_parse
    implicit none
    private

    public :: parse_toml_info

contains

    !> 解析程序所需的所有配置文件
    subroutine parse_toml_info()
        call parse_access_toml()
        call parse_info_sph()
    end subroutine parse_toml_info

    !> 解析access.toml文件
    !>
    !> 用法: fpm run sph -- example/shearcavity
    subroutine parse_access_toml()
        type(toml_table), allocatable :: access_table
        type(toml_table), pointer :: subtable
        character(:), allocatable :: path
        integer access_toml_unit

        call get_path(path)
        if (.not. is_exist(path//"/access.toml")) &
            call error_stop(path//'/access.toml文件不存在, 请检查!', &
                            'toml_info_m%parse_access_toml')
        open (newunit=access_toml_unit, file=path//'/access.toml', status='old')
        call toml_parse(access_table, access_toml_unit)
        close (access_toml_unit)

        call get_value(access_table, 'access', subtable)
        call get_value(subtable, 'in', in_path, './data/input')  !! 默认输入路径
        call get_value(subtable, 'out', out_path, './data/output')  !! 默认输出路径

        ! 日志
        call stdlog%log_information('In  path: '//in_path)
        call stdlog%log_information('Out path: '//out_path)

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
        call get_value(subtable, 'virtual_part', virtual_part, .true.)
        call get_value(subtable, 'self_gravity', self_gravity, .false.)
        call get_value(subtable, 'visc', visc, .true.)
        call get_value(subtable, 'eos_form', eos_form, 2)
        call get_value(subtable, 'visc_artificial', visc_artificial, .false.)
        call get_value(subtable, 'heat_artificial', heat_artificial, .false.)

        call get_value(sph_table, 'pre-process', subtable)
        call get_value(subtable, 'dofile', dofile, .false.) ! 默认不从 Lua 脚本中生成数据
        call get_value(subtable, 'lua_script', lua_script, 'lua_script.lua')

        ! 日志
        call stdlog%log_information('Project name: '//nick)
        call stdlog%log_information('dt: '//to_string(dt))
        call stdlog%log_information('skf: '//to_string(skf))
        call stdlog%log_information('nnps: '//to_string(nnps))
        call stdlog%log_information('print_step: '//to_string(print_step))
        call stdlog%log_information('save_step: '//to_string(save_step))
        call stdlog%log_information('maxn: '//to_string(maxn))
        call stdlog%log_information('kpair: '//to_string(kpair))
        call stdlog%log_information('virtual_part: '//to_string(virtual_part))
        call stdlog%log_information('self_gravity: '//to_string(self_gravity))
        call stdlog%log_information('visc: '//to_string(visc))
        call stdlog%log_information('eos_form: '//to_string(eos_form))
        call stdlog%log_information('visc_artificial: '//to_string(visc_artificial))
        call stdlog%log_information('heat_artificial: '//to_string(heat_artificial))
        call stdlog%log_information('dofile: '//to_string(dofile))
        call stdlog%log_information('lua_script: '//lua_script)

        nullify (subtable)
    end subroutine parse_info_sph

end module toml_info_m
