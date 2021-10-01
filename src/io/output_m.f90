module output_m

    use iso_fortran_env, only: int8
    use utils, only: to_string
    implicit none

contains

    !> 输出每个保存时间步的求解信息（拓展）
    !>     subroutine for saving particle information to external disk file
    !>
    !>     x-- coordinates of particles                                  [in]
    !>     vx-- velocities of particles                                  [in]
    !>     mass-- mass of particles                                      [in]
    !>     rho-- dnesities of particles                                  [in]
    !>     p-- pressure  of particles                                    [in]
    !>     u-- internal energy of particles                              [in]
    !>     c-- sound velocity of particles                               [in]
    !>     itype-- types of particles                                    [in]
    !>     hsml-- smoothing lengths of particles                         [in]
    !>     ntotal-- total particle number                                [in]

    subroutine output_all(x, vx, mass, rho, p, u, c, itype, hsml, ntotal, n)

        use sph_kind, only: rk
        use parameter
        implicit none

        real(rk), intent(in) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), u(maxn), c(maxn), hsml(maxn)
        integer, intent(in)  :: itype(maxn), ntotal, n
        
        integer  :: i, d, npart
        integer(int8) :: xv_unit, state_unit, other_unit

        !> 输出粒子的位置、速度信息
        open (newunit=xv_unit   , file='./example/data/all/f_'//to_string(n)//'xv.dat'   )
        
        !> 输出粒子的宏观信息：质量、密度、压强、内能
        open (newunit=state_unit, file='./example/data/all/f_'//to_string(n)//'state.dat')

        !> 输出粒子的其它信息：粒子类型、光滑长度
        open (newunit=other_unit, file='./example/data/all/f_'//to_string(n)//'other.dat')

        write (xv_unit, *) ntotal
        do i = 1, ntotal
            write (xv_unit   , 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            write (state_unit, 1002) i, mass(i), rho(i), p(i), u(i)
            write (other_unit, 1003) i, itype(i), hsml(i)
        end do

        close (xv_unit   )
        close (state_unit)
        close (other_unit)

    1001 format(1x, i6, 6(2x, e14.8))
    1002 format(1x, i6, 7(2x, e14.8))
    1003 format(1x, i6, 2x, i4, 2x, e14.8)

    end subroutine output_all


end module output_m