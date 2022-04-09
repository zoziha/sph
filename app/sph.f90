!> 主程序，读入配置文件，初始化数据，启动计算服务。
!>
!> - 作者: 左志华
!> - 日期: 2022年3月
program main

    use config_m, only: rk, stdout, stdin, tinsert, tsearch, dt, nnps, maxn, dofile, &
        lua_script
    use info_m, only: operator(.c.), info
    use input_m, only: input
    use lua_call_m, only: lua_input
    use macro_m, only: x, vx, mass, rho, p, u, c, s, e, hsml, itype, alloc_macro_memory
    use master_time_m, only: tic, toc, time_print
    use output_m, only: set_parameter_log, set_folder
    use parameter, only: dim
    use stdlib_logger, only: stdlog => global_logger
    use time_integration_m, only: time_integration
    use toml_info_m, only: parse_toml_info
    implicit none

    integer :: ntotal       !! 粒子总数
    integer :: maxtimestep  !! 最大时间步数
    integer :: yesorno      !! 是否继续运行

    call tic()
    call time_print()
    call stdlog%add_log_file('.stdlog.log')
    call stdlog%log_information('Start logging')

    ! 前处理
    call parse_toml_info()
    call alloc_macro_memory(dim, maxn, x, vx, mass, rho, p, u, c, s, e, hsml, itype)
    call set_parameter_log()
    call set_folder()

    if (dofile) then
        call lua_input(lua_script, x, vx, mass, rho, p, u, itype, hsml, ntotal)
        call stdlog%log_information('Loaded lua script successfully')
    else
        call input(x, vx, mass, rho, p, u, itype, hsml, ntotal)
        call stdlog%log_information('Loaded non-lua input successfully')
    end if

    ! 主循环
    do
        write (stdout, "(a)", advance="no") .c.'<c>Please input the maximal time steps: '
        read (stdin, *) maxtimestep

        if (maxtimestep > 0) &
            call time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt)

        write (stdout, "(a)", advance="no") .c.'<c>Are you going to run more time steps ? (0=no, 1=yes): '
        read (stdin, *) yesorno
        if (yesorno == 0) exit

    end do

    if (nnps == 3) write (stdout, '(2(a,es10.3),a)') info('Particle insertion time: '), tinsert, &
        's, particle search time: ', tsearch, 's'

    call stdlog%log_information('End logging')
    call time_print()
    call toc()
    write (stdout, "(a)") 'All finish!'

end program main
