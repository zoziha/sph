!> ParaView 辅助程序
!> 属于 SPH 下属项目
!>
!> - 作者: 左志华
!> - 日期: 2022-3-31
program to_vtk

    use paraview_interface, only: output_to_paraview_vtk
    implicit none

    call output_to_paraview_vtk()

end program to_vtk
