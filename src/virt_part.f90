!>   subroutine to determine the information of virtual particles
!>   here only the monaghan type virtual particles for the 2d shear
!>   cavity driven problem are generated.
!>
!>     itimestep : current time step                                 [in]
!>     ntotal : number of particles                                  [in]
!>     nvirt  : number of virtual particles                         [out]
!>     hsml   : smoothing length                                 [in|out]
!>     mass   : particle masses                                  [in|out]
!>     x      : coordinates of all particles                     [in|out]
!>     vx     : velocities of all particles                      [in|out]
!>     rho    : density                                          [in|out]
!>     u      : internal energy                                  [in|out]
!>     itype   : type of particles                               [in|out]

subroutine virt_part(itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype)

    use sph_kind, only: rk
    use parameter
    implicit none

    integer  :: itimestep, ntotal, nvirt, itype(maxn)
    real(rk) :: hsml(maxn), mass(maxn), x(dim, maxn), vx(dim, maxn), rho(maxn), u(maxn), p(maxn)
    integer  :: i, j, d, im, mp
    real(rk) :: xl, dx, v_inf

    if (vp_input) then

        open (1, file='./example/data/xv_vp.dat')
        open (2, file='./example/data/state_vp.dat')
        open (3, file='./example/data/other_vp.dat')
        read (1, *) nvirt
        do j = 1, nvirt
            i = ntotal + j
            read (1, *) im, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            read (2, *) im, mass(i), rho(i), p(i), u(i)
            read (3, *) im, itype(i), hsml(i)
        end do
        close (1)
        close (2)
        close (3)

    else

        nvirt = 0
        mp    = 40
        xl    = 1.0e-3
        dx    = xl/mp
        v_inf = 1.e-3

        !     monaghan type virtual particle on the upper side

        do i = 1, 2*mp + 1
            nvirt = nvirt + 1
            x(1, ntotal + nvirt)  = (i - 1)*dx/2
            x(2, ntotal + nvirt)  = xl
            vx(1, ntotal + nvirt) = v_inf
            vx(2, ntotal + nvirt) = 0.
        end do

        !     monaghan type virtual particle on the lower side

        do i = 1, 2*mp + 1
            nvirt = nvirt + 1
            x(1, ntotal + nvirt)  = (i - 1)*dx/2
            x(2, ntotal + nvirt)  = 0.
            vx(1, ntotal + nvirt) = 0.
            vx(2, ntotal + nvirt) = 0.
        end do

        !     monaghan type virtual particle on the left side

        do i = 1, 2*mp - 1
            nvirt = nvirt + 1
            x(1, ntotal + nvirt)  = 0.
            x(2, ntotal + nvirt)  = i*dx/2
            vx(1, ntotal + nvirt) = 0.
            vx(2, ntotal + nvirt) = 0.
        end do

        !     monaghan type virtual particle on the right side

        do i = 1, 2*mp - 1
            nvirt = nvirt + 1
            x(1, ntotal + nvirt)  = xl
            x(2, ntotal + nvirt)  = i*dx/2
            vx(1, ntotal + nvirt) = 0.
            vx(2, ntotal + nvirt) = 0.
        end do

        do i = 1, nvirt
            rho(ntotal + i)   = 1000.
            mass(ntotal + i)  = rho(ntotal + i)*dx*dx
            p(ntotal + i)     = 0.
            u(ntotal + i)     = 357.1
            itype(ntotal + i) = -2
            hsml(ntotal + i)  = dx
        end do

    end if

    if (mod(itimestep, save_step) == 0) then
        open (1, file='./example/data/xv_vp.dat')
        open (2, file='./example/data/state_vp.dat')
        open (3, file='./example/data/other_vp.dat')
        write (1, *) nvirt
        do i = ntotal + 1, ntotal + nvirt
            write (1, 1001) i, (x(d, i), d=1, dim), (vx(d, i), d=1, dim)
            write (2, 1002) i, mass(i), rho(i), p(i), u(i)
            write (3, 1003) i, itype(i), hsml(i)
        end do
        close (1)
        close (2)
        close (3)
    end if

    if (mod(itimestep, print_step) == 0) then
        if (int_stat) then
            print *, ' >> statistics: virtual boundary particles:'
            print *, '          number of virtual particles:', nvirt
        end if
    end if
    
1001 format(1x, i6, 6(2x, e14.8))
1002 format(1x, i6, 7(2x, e14.8))
1003 format(1x, i6, 2x, i4, 2x, e14.8)

end subroutine virt_part
