!>     subroutine for loading or generating initial particle information
!>
!>     x-- coordinates of particles                                 [out]
!>     vx-- velocities of particles                                 [out]
!>     mass-- mass of particles                                     [out]
!>     rho-- dnesities of particles                                 [out]
!>     p-- pressure  of particles                                   [out]
!>     u-- internal energy of particles                             [out]
!>     itype-- types of particles                                   [out]
!>     hsml-- smoothing lengths of particles                        [out]
!>     ntotal-- total particle number                               [out]

subroutine input(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    use sph_kind, only: rk
    implicit none
    include 'param.inc'

    integer  :: itype(maxn), ntotal
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), p(maxn), u(maxn), hsml(maxn), rho(maxn)
    integer  :: i, d, im

    !     load initial particle information from external disk file

    if (config_input) then

        open (1, file='./example/data/f_xv.dat')
        open (2, file='./example/data/f_state.dat')
        open (3, file='./example/data/f_other.dat')

        write (*, *) '  **************************************************'
        write (*, *) '      loading initial particle configuration...   '
        read (1, *) ntotal
        write (*, *) '      total number of particles   ', ntotal
        write (*, *) '  **************************************************'
        do i = 1, ntotal
            read (1, *) im, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            read (2, *) im, mass(i), rho(i), p(i), u(i)
            read (3, *) im, itype(i), hsml(i)
        end do

    else

        open (1, file='./example/data/ini_xv.dat')
        open (2, file='./example/data/ini_state.dat')
        open (3, file='./example/data/ini_other.dat')

        if (shocktube) call shock_tube(x, vx, mass, rho, p, u, itype, hsml, ntotal)

        if (shearcavity) call shear_cavity(x, vx, mass, rho, p, u, itype, hsml, ntotal)
        do i = 1, ntotal
            write (1, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            write (2, 1002) i, mass(i), rho(i), p(i), u(i)
            write (3, 1003) i, itype(i), hsml(i)
        end do
        write (*, *) '  **************************************************'
        write (*, *) '      initial particle configuration generated   '
        write (*, *) '      total number of particles   ', ntotal
        write (*, *) '  **************************************************'

    end if

    close (1)
    close (2)
    close (3)

1001 format(1x, i5, 6(2x, e14.8))
1002 format(1x, i5, 7(2x, e14.8))
1003 format(1x, i5, 2x, i2, 2x, e14.8)

end subroutine input

!>     this subroutine is used to generate initial data for the
!>     1 d noh shock tube problem
!>     x-- coordinates of particles                                 [out]
!>     vx-- velocities of particles                                 [out]
!>     mass-- mass of particles                                     [out]
!>     rho-- dnesities of particles                                 [out]
!>     p-- pressure  of particles                                   [out]
!>     u-- internal energy of particles                             [out]
!>     itype-- types of particles                                   [out]
!>          =1   ideal gas
!>     hsml-- smoothing lengths of particles                        [out]
!>     ntotal-- total particle number                               [out]

subroutine shock_tube(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    use sph_kind, only: rk
    implicit none
    include 'param.inc'

    integer  :: itype(maxn), ntotal
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), u(maxn), hsml(maxn)
    integer  :: i, d
    real(rk) :: space_x

    ntotal = 400
    space_x = 0.6/80.

    do i = 1, ntotal
        mass(i) = 0.75/400.
        hsml(i) = 0.015
        itype(i) = 1
        do d = 1, dim
            x(d, i) = 0.
            vx(d, i) = 0.
        end do
    end do

    do i = 1, 320
        x(1, i) = -0.6 + space_x/4.*(i - 1)
    end do

    do i = 320 + 1, ntotal
        x(1, i) = 0.+space_x*(i - 320)
    end do

    do i = 1, ntotal
        if (x(1, i) <= 1.e-8) then
            u(i) = 2.5
            rho(i) = 1.
            p(i) = 1.
        end if
        if (x(1, i) > 1.e-8) then
            u(i) = 1.795
            rho(i) = 0.25
            p(i) = 0.1795
        end if
    end do

end subroutine shock_tube

subroutine shear_cavity(x, vx, mass, rho, p, u, itype, hsml, ntotal)

    !----------------------------------------------------------------------
    !     this subroutine is used to generate initial data for the
    !     2 d shear driven cavity probem with re = 1
    !     x-- coordinates of particles                                 [out]
    !     vx-- velocities of particles                                 [out]
    !     mass-- mass of particles                                     [out]
    !     rho-- dnesities of particles                                 [out]
    !     p-- pressure  of particles                                   [out]
    !     u-- internal energy of particles                             [out]
    !     itype-- types of particles                                   [out]
    !          =2   water
    !     h-- smoothing lengths of particles                           [out]
    !     ntotal-- total particle number                               [out]

    use sph_kind, only: rk
    implicit none
    include 'param.inc'

    integer  :: itype(maxn), ntotal
    real(rk) :: x(dim, maxn), vx(dim, maxn), mass(maxn), rho(maxn), p(maxn), u(maxn), hsml(maxn)
    integer  :: i, j, d, m, n, mp, np, k
    real(rk) :: xl, yl, dx, dy

    !     giving mass and smoothing length as well as other data.

    m = 41
    n = 41
    mp = m - 1
    np = n - 1
    ntotal = mp*np
    xl = 1.e-3
    yl = 1.e-3
    dx = xl/mp
    dy = yl/np

    do i = 1, mp
        do j = 1, np
            k = j + (i - 1)*np
            x(1, k) = (i - 1)*dx + dx/2.
            x(2, k) = (j - 1)*dy + dy/2.
        end do
    end do

    do i = 1, mp*np
        vx(1, i) = 0.
        vx(2, i) = 0.
        rho(i) = 1000.
        mass(i) = dx*dy*rho(i)
        p(i) = 0.
        u(i) = 357.1
        itype(i) = 2
        hsml(i) = dx
    end do

end subroutine shear_cavity
