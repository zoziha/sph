!> 将数据导入paraview(简易接口)
module paraview_interface

    use parameter
    use utils, only: to_string
    implicit none
    private

    public :: output_to_paraview_vtk

contains

    !> 将数据输出为paraview可读的vtk格式
    subroutine output_to_paraview_vtk()
        integer :: steps, i

        print *, "请输入时间步数："
        read (*, *) steps

        do i = 1, steps
            print "(a, i0, a)", "正在输出第", i, "个时间步。"
            call output_to_paraview_vtk_one_step(i)
        end do

    end subroutine output_to_paraview_vtk

    subroutine output_to_paraview_vtk_one_step(i_steps)
        integer, intent(in) :: i_steps
        integer :: ntotal
        integer :: unit_xv, unit_state, unit_other
        integer :: index, i, d
        real(rk), allocatable :: x(:, :), vx(:, :), mass(:), rho(:), p(:), u(:), itype(:), hsml(:)

        !> 读入数据:xv, state, other
        open (newunit=unit_xv, file="example/data/all/f_"//to_string(i_steps)//"xv.dat", status="old")
        open (newunit=unit_state, file="example/data/all/f_"//to_string(i_steps)//"state.dat", status="old")
        open (newunit=unit_other, file="example/data/all/f_"//to_string(i_steps)//"other.dat", status="old")

        read (unit_xv, *) ntotal
        allocate (x(dim, ntotal), vx(dim, ntotal), mass(ntotal), rho(ntotal), p(ntotal), u(ntotal), itype(ntotal), hsml(ntotal))

        do i = 1, ntotal
            read (unit_xv, *) index, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            read (unit_state, *) index, mass(i), rho(i), p(i), u(i)
            read (unit_other, *) index, itype(i), hsml(i)
        end do

        close (unit_xv)
        close (unit_state)
        close (unit_other)

        !> 输出数据: vtk
        open (newunit=unit_xv, file="example/data/paraview/sph"//to_string(i_steps)//".vtk", access="stream", &
              form="formatted")

        !> 输出头部和点坐标
        write (unit_xv, "(a)") "# vtk DataFile Version 3.0"
        write (unit_xv, "(a)") "paraview_vtk_output"
        write (unit_xv, "(a)") "ASCII"
        write (unit_xv, "(a)") "DATASET UNSTRUCTURED_GRID"
        write (unit_xv, '(a, i0, a)') "POINTS ", ntotal, " float"
        do i = 1, ntotal
            if (dim == 2) then
                write (unit_xv, "(*(ES12.5, 3x))") x(:, i), 0.0  !! 第三个坐标不能缺少
            elseif (dim == 3) then
                write (unit_xv, "(*(ES12.5, 3x))") x(:, i)
            end if
        end do
        write (unit_xv, "(a, i0)") "POINT_DATA ", ntotal

        !> 输出点的质量属性
        write (unit_xv, "(a)") "SCALARS mass float 1"
        write (unit_xv, "(a)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (unit_xv, "(*(ES12.5, 3x))") mass(i)
        end do

        !> 输出点的密度属性
        write (unit_xv, "(a)") "SCALARS rho float 1"
        write (unit_xv, "(a)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (unit_xv, "(*(ES12.5, 3x))") rho(i)
        end do

        !> 输出点的压力属性
        write (unit_xv, "(a)") "SCALARS p float 1"
        write (unit_xv, "(a)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (unit_xv, "(*(ES12.5, 3x))") p(i)
        end do

        !> 输出点的内能属性
        write (unit_xv, "(a)") "SCALARS u float 1"
        write (unit_xv, "(a)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (unit_xv, "(*(ES12.5, 3x))") u(i)
        end do

        !> 输出点的类型属性
        write (unit_xv, "(a)") "SCALARS itype float 1"
        write (unit_xv, "(a)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (unit_xv, "(ES12.5, 3x)") itype(i)
        end do

        !> 输出点的光滑长度属性
        write (unit_xv, "(a)") "SCALARS hsml float 1"
        write (unit_xv, "(a)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (unit_xv, "(*(ES12.5, 3x))") hsml(i)
        end do

        !> 输出点的速度属性
        write (unit_xv, "(a)") "VECTORS vx float"
        do i = 1, ntotal
            if (dim == 2) then
                write (unit_xv, "(*(ES12.5, 3x))") vx(:, i), 0.0
            elseif (dim == 3) then
                write (unit_xv, "(*(ES12.5, 3x))") vx(:, i)
            end if
        end do
        close (unit_xv)

    end subroutine output_to_paraview_vtk_one_step

end module paraview_interface
