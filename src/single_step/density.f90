!> 通过应用密度求和法更新密度的子程序。详见第 4 章中的论述。
!>   subroutine to calculate the density with sph summation algorithm.
!>   see equ.(4.35)
!>
!>     ntotal : number of particles                                  [in]
!>     hsml   : smoothing length                                     [in]
!>     mass   : particle masses                                      [in]
!>     niac   : number of interaction pairs                          [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     w      : kernel for all interaction pairs                     [in]
!>     itype   : type of particles                                   [in]
!>     x       : coordinates of all particles                        [in]
!>     rho    : density                                             [out]

subroutine sum_density(ntotal, hsml, mass, niac, pair_i, pair_j, w, itype, rho)

    use sph_kinds, only: rk
    use parameter
    implicit none

    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 相互作用对的数目
    integer, intent(in) :: niac
    integer :: pair_i(max_interaction), pair_j(max_interaction), itype(maxn)
    real(rk) :: hsml(maxn), mass(maxn), w(max_interaction), rho(maxn)
    integer :: i, j, k
    real(rk) :: selfdens, hv(dim), wi(maxn)

    !     wi(maxn)---integration of the kernel itself

    hv(1:dim) = 0.0_rk

    ! 先计算粒子自身对密度的贡献值
    !     self density of each particle: wii (kernel for distance 0)
    !     and take contribution of particle itself:

    ! 修正密度加和法, 先计算分母
    !     firstly calculate the integration of the kernel over the space

    if (nor_density) then
        do i = 1, ntotal
            call kernel(0.0_rk, hv, hsml(i), selfdens, hv)
            wi(i) = selfdens*mass(i)/rho(i)
        end do

        do k = 1, niac
            i = pair_i(k)
            j = pair_j(k)
            wi(i) = wi(i) + mass(j)/rho(j)*w(k)
            wi(j) = wi(j) + mass(i)/rho(i)*w(k)
        end do

    end if
    ! 修正密度加和法, 再计算分子
    !     secondly calculate the rho integration over the space

    do i = 1, ntotal
        call kernel(0.0_rk, hv, hsml(i), selfdens, hv)
        rho(i) = selfdens*mass(i)
    end do

    !     calculate sph sum for rho:
    do k = 1, niac
        i = pair_i(k)
        j = pair_j(k)
        rho(i) = rho(i) + mass(j)*w(k)
        rho(j) = rho(j) + mass(i)*w(k)
    end do

    ! 修正密度加和法, 正则化
    !     thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

    if (nor_density) then
        rho(i:ntotal) = rho(i:ntotal)/wi(i:ntotal)
    end if

end subroutine sum_density

!> 通过应用连续密度法更新密度的子程序。详见第 4 章中的论述。
!>     subroutine to calculate the density with sph continuity approach.
!>     see equ.(4.34)
!>
!>     ntotal : number of particles                                  [in]
!>     mass   : particle masses                                      [in]
!>     niac   : number of interaction pairs                          [in]
!>     pair_i : list of first partner of interaction pair            [in]
!>     pair_j : list of second partner of interaction pair           [in]
!>     dwdx   : derivation of kernel for all interaction pairs       [in]
!>     vx     : velocities of all particles                          [in]
!>     itype   : type of particles                                   [in]
!>     x      : coordinates of all particles                         [in]
!>     rho    : density                                              [in]
!>     drhodt : density change rate of each particle                [out]

subroutine con_density(ntotal, mass, niac, pair_i, pair_j, dwdx, vx, itype, x, rho, drhodt)

    use sph_kinds, only: rk
    use parameter
    implicit none

    !> 在模拟中所使用的粒子总数
    !> number of particles in simulation
    integer, intent(in) :: ntotal
    !> 相互作用对的数目
    integer, intent(in) :: niac
    integer :: pair_i(max_interaction), pair_j(max_interaction), itype(maxn)
    real(rk) :: mass(maxn), dwdx(dim, max_interaction), vx(dim, maxn), x(dim, maxn), rho(maxn), drhodt(maxn)
    integer :: i, j, k, d
    real(rk) :: vcc, dvx(dim)

    do i = 1, ntotal
        drhodt(i) = 0._rk
    end do

    do k = 1, niac

        i = pair_i(k)
        j = pair_j(k)
        do d = 1, dim
            dvx(d) = vx(d, i) - vx(d, j)
        end do
        vcc = dvx(1)*dwdx(1, k)
        do d = 2, dim
            vcc = vcc + dvx(d)*dwdx(d, k)
        end do
        drhodt(i) = drhodt(i) + mass(j)*vcc
        drhodt(j) = drhodt(j) + mass(i)*vcc

    end do

end subroutine con_density
