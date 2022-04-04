! 主程序，读入配置文件，初始化数据，启动计算服务。
! -- 作者: 左志华
! -- 日期: 2022年3月
program main

    use config_m, only: rk, stdout, stdin, tinsert, tsearch, dt, nnps
    use parameter
    use master_time_m, only: tic, toc, time_print
    use output_m, only: set_parameter_log, set_folder
    use toml_info_m, only: parse_toml_info
    use macro_m, only: x, vx, mass, rho, p, u, c, s, e, hsml, itype, alloc_macro_memory
    use info_m, only: operator(.c.), info
    use input_m, only: input
    use time_integration_m, only: time_integration
    implicit none

    integer :: ntotal
    integer :: maxtimestep
    integer :: yesorno

    call tic()
    call time_print()

    ! 前处理
    call parse_toml_info()
    call alloc_macro_memory(dim, maxn, x, vx, mass, rho, p, u, c, s, e, hsml, itype)
    call set_parameter_log()
    call set_folder()

    call input(x, vx, mass, rho, p, u, itype, hsml, ntotal)

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

    call time_print()
    call toc()
    write (stdout, "(a)") 'All finish!'

end program main
