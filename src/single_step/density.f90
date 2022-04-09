!> 密度
module density_m

    use config_m, only: rk
    use parameter
    use kernel_m, only: kernel
    implicit none
    private

    public :: sum_density, con_density
    
contains

    !> 通过应用密度求和法更新密度的子程序。详见第 4 章中的论述 (式4.35)。
    subroutine sum_density(ntotal, hsml, mass, niac, pair_i, pair_j, w, itype, rho)
        integer, intent(in) :: ntotal       !! 在模拟中所使用的粒子总数
        real(rk), intent(in) :: hsml(:)     !! 光滑长度
        real(rk), intent(in) :: mass(:)     !! 粒子质量
        integer, intent(in) :: niac         !! 相互作用对的数目
        integer, intent(in) :: pair_i(:)    !! 相互作用对的第一个粒子
        integer, intent(in) :: pair_j(:)    !! 相互作用对的第二个粒子
        real(rk), intent(in) :: w(:)        !! 相互作用对的核函数
        integer, intent(in) :: itype(:)     !! 粒子类型
        real(rk), intent(out) :: rho(:)     !! 密度
        
        integer :: i, j, k
        real(rk) :: selfdens, hv(dim), wi(ntotal)

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

    !> 通过应用连续密度法更新密度的子程序。详见第 4 章中的论述 (式4.34)
    pure subroutine con_density(ntotal, mass, niac, pair_i, pair_j, dwdx, vx, itype, x, rho, drhodt)
        integer, intent(in) :: ntotal       !! 在模拟中所使用的粒子总数
        real(rk), intent(in) :: mass(:)     !! 粒子质量
        integer, intent(in) :: niac         !! 相互作用对的数目
        integer, intent(in) :: pair_i(:)    !! 相互作用对的第一个粒子
        integer, intent(in) :: pair_j(:)    !! 相互作用对的第二个粒子
        real(rk), intent(in) :: dwdx(:, :)  !! 相互作用对的核函数
        real(rk), intent(in) :: vx(:, :)    !! 粒子速度
        integer, intent(in) :: itype(:)     !! 粒子类型
        real(rk), intent(in) :: x(:, :)     !! 粒子坐标
        real(rk), intent(in) :: rho(:)      !! 密度
        real(rk), intent(out) :: drhodt(:)  !! 密度变化率
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
    
end module density_m
