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

subroutine output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)



    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: itype(maxn), ntotal
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), u(maxn), c(maxn), hsml(maxn)
    integer  :: i, d, npart

    open (1, file='./example/data/f_xv.dat')
    open (2, file='./example/data/f_state.dat')
    open (3, file='./example/data/f_other.dat')

    write (1, *) ntotal
    do i = 1, ntotal
        write (1, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
        write (2, 1002) i, mass(i), rho(i), p(i), u(i)
        write (3, 1003) i, itype(i), hsml(i)
    end do

    close (1)
    close (2)
    close (3)

1001 format(1x, i6, 6(2x, e14.8))
1002 format(1x, i6, 7(2x, e14.8))
1003 format(1x, i6, 2x, i4, 2x, e14.8)

end subroutine output
