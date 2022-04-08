!> 输出模块
module output_m

    use easy_string_m, only: to_string
    use swift_file_m, only: mkdir, is_exist
    use, intrinsic :: iso_fortran_env, only: stdout => output_unit
    use config_m, only: nick, out_path, skf, nnps, print_step
    use info_m, only: operator(.c.)
    use error_stop_m, only: error_stop
    use parameter
    implicit none

    public :: set_folder, set_parameter_log

contains

    !> 输出每个保存时间步的求解信息（拓展）
    !> subroutine for saving particle information to external disk file
    subroutine output_all(x, vx, mass, rho, p, u, c, itype, hsml, ntotal, n)

        use config_m, only: rk
        use parameter
        implicit none

        !> 粒子的位置
        !> coordinates of particles
        real(rk), intent(in) :: x(:, :)
        !> 粒子的速度
        !> velocities of particles
        real(rk), intent(in) :: vx(:, :)
        !> 粒子的质量
        !> mass of particles
        real(rk), intent(in) :: mass(:)
        !> 粒子的密度
        !> dnesities of particles
        real(rk), intent(in) :: rho(:)
        !> 粒子的压力
        !> pressure  of particles
        real(rk), intent(in) :: p(:)
        !> 粒子的内部能量
        !> internal energy of particles
        real(rk), intent(in) :: u(:)
        !> 粒子的声速
        !> sound velocity of particles
        real(rk), intent(in) :: c(:)
        !> 粒子的类型(1: ideal gas; 2: water; 3: TNT)
        !> types of particles
        integer, intent(in) :: itype(:)
        !> 粒子的光滑长度
        !> smoothing lengths of particles
        real(rk), intent(in) :: hsml(:)
        !> 在模拟中所使用的粒子总数
        !> number of particles in simulation
        integer, intent(in) :: ntotal
        !> 第几个时间步
        !> number of time step
        integer, intent(in) :: n

        integer :: i, d
        integer :: xv_unit, state_unit, other_unit

        !> 输出粒子的位置、速度信息
        open (newunit=xv_unit, file=out_path//'/all/f_'//to_string(n)//'xv.dat')

        !> 输出粒子的宏观信息：质量、密度、压强、内能
        open (newunit=state_unit, file=out_path//'/all/f_'//to_string(n)//'state.dat')

        !> 输出粒子的其它信息：粒子类型、光滑长度
        open (newunit=other_unit, file=out_path//'/all/f_'//to_string(n)//'other.dat')

        write (xv_unit, *) ntotal
        do i = 1, ntotal
            write (xv_unit, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            write (state_unit, 1002) i, mass(i), rho(i), p(i), u(i)
            write (other_unit, 1003) i, itype(i), hsml(i)
        end do

        close (xv_unit)
        close (state_unit)
        close (other_unit)

1001    format(1x, i6, 6(2x, e14.8))
1002    format(1x, i6, 7(2x, e14.8))
1003    format(1x, i6, 2x, i4, 2x, e14.8)

    end subroutine output_all

    !> 输出粒子作用对的统计信息
    subroutine set_statistics_print(itimestep, ntotal, niac, countiac)
        integer, intent(in) :: itimestep, ntotal, niac, countiac(:)
        integer sumiac, maxiac, miniac, noiac, maxp, minp, i

        sumiac = 0
        maxiac = 0
        miniac = 1000
        noiac = 0
        do concurrent(i=1:ntotal)
            sumiac = sumiac + countiac(i)
            if (countiac(i) > maxiac) then
                maxiac = countiac(i)
                maxp = i
            end if
            if (countiac(i) < miniac) then
                miniac = countiac(i)
                minp = i
            end if
            if (countiac(i) == 0) noiac = noiac + 1
        end do

        if (mod(itimestep, print_step) == 0) then
            if (int_stat) then
                write (stdout, '(/a)') .c.'Statistics: interactions per particle:'
                print 100, 'particle: ', maxp, ' maximal interactions: ', maxiac
                print 100, 'particle: ', minp, ' minimal interactions: ', miniac
                ! 平均每个粒子的作用对数
                print 100, 'average : ', sumiac/ntotal
                print 100, 'total pairs : ', niac
                print 100, 'particles with no interactions: ', noiac
            end if
        end if

100     format(*(a, i0))
    end subroutine set_statistics_print

    !> 输出工程特征信息
    subroutine set_parameter_log()
        write (stdout, '(a)') .c.'Project name: '//nick
        write (stdout, '(a,i0)') .c.'Smoothed kernel function: ', skf
        write (stdout, '(a,i0/)') .c.'NNPS method: ', nnps
    end subroutine set_parameter_log

    !> 建立所需文件夹
    subroutine set_folder()
        !@tofix: mkdir 只能生成一级目录，不能生成多级目录
        if (.not. is_exist(out_path, .true.)) call mkdir(out_path, .true.)
        if (.not. is_exist(out_path, .true.)) call error_stop('cannot create folder: '//out_path, &
                                                              'output_m%set_folder')
        if (.not. is_exist(out_path//'/all', .true.)) call mkdir(out_path//'/all', .true.)
        if (.not. is_exist(out_path//'/paraview', .true.)) call mkdir(out_path//'/paraview', .true.)
        if (.not. is_exist(out_path//'/all', .true.)) call mkdir(out_path//'/all', .true.)
    end subroutine set_folder

end module output_m
