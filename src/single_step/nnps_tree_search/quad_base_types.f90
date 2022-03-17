module quad_types

    use config_m, only: rk
    implicit none
    private
    
    public :: point_t, rectangle_t, circle_t

    !> 二维点
    type point_t
        real(rk) :: x  !! 点的横坐标
        real(rk) :: y  !! 点的纵坐标
        integer :: index !! 粒子索引
    end type point_t

    !> 矩形
    type rectangle_t

        real(rk) :: x  !! 矩形形心横坐标
        real(rk) :: y  !! 矩形形心纵坐标
        real(rk) :: w  !! 矩形宽度
        real(rk) :: h  !! 矩形高度

    contains

        procedure :: contains => rectangle_t_contains
        procedure :: intersects => rectangle_t_intersects

    end type rectangle_t

    !> 圆形
    type circle_t

        real(rk) :: x, y, r

    contains

        procedure :: contains => circle_t_contains
        procedure :: intersects => circle_t_intersects

    end type circle_t

contains

    !> 查询粒子是否在矩形域内
    !> 细分算法会导致某些在矩形边上的粒子无法插入树型表
    logical function rectangle_t_contains(self, point) result(contains)

        class(rectangle_t), intent(in) :: self
        type(point_t), intent(in) :: point

        associate (x => self%x, &
                   y => self%y, &
                   w => self%w, &
                   h => self%h)

            contains = (point%x > x - 0.5_rk*w) .and. &
                       (point%x < x + 0.5_rk*w) .and. &
                       (point%y > y - 0.5_rk*h) .and. &
                       (point%y < y + 0.5_rk*h)

        end associate

    end function rectangle_t_contains

    !> 查询几何形状是否有交集
    logical function rectangle_t_intersects(self, range) result(intersects)

        class(rectangle_t), intent(in) :: self
        type(rectangle_t), intent(in) :: range

        !&< fprettify flag
        associate (left    => self%x  - 0.5_rk*self%w , &
                   right   => self%x  + 0.5_rk*self%w , &
                   bottom  => self%y  - 0.5_rk*self%h , &
                   top     => self%y  + 0.5_rk*self%h , &
                   left_   => range%x - 0.5_rk*range%w, &
                   right_  => range%x + 0.5_rk*range%w, &
                   bottom_ => range%y - 0.5_rk*range%h, &
                   top_    => range%y + 0.5_rk*range%h)

            intersects = (left   < right_) .or. &
                         (right  > left_ ) .or. &
                         (bottom < top_  ) .or. &
                         (top    > bottom_)

        end associate
        !&>

    end function rectangle_t_intersects

    !> 查询粒子是否在圆形内
    logical function circle_t_contains(self, point) result(contains)

        class(circle_t), intent(in) :: self
        type(point_t), intent(in) :: point

        contains = hypot(self%x - point%x, self%y - point%y) < self%r

    end function circle_t_contains

    !> 查询几何形状是否有交集
    logical function circle_t_intersects(self, range) result(intersects)

        class(circle_t), intent(in) :: self
        type(rectangle_t), intent(in) :: range

        real(rk) :: x_dist, y_dist

        ! 形心距离
        x_dist = abs(range%x - self%x)
        y_dist = abs(range%y - self%y)

        !&< fprettify flag
        associate (r  => self%r     , &
                   w_ => 0.5*range%w, &
                   h_ => 0.5*range%h)
            associate (edges => hypot(x_dist - w_, y_dist - h_))

                !> 圆形边缘在矩形外
                !> no intersection
                if (x_dist >= r + w_ .or. y_dist >= r + h_) then
                    intersects = .false.
                    return
                end if

                !> 圆形形心在矩形内
                !> intersection within the circle
                if (x_dist < w_ .or. y_dist < h_) then
                    intersects = .true.
                    return
                end if

                !> 存在矩形顶点在圆形内
                !> intersection on the edge of the circle
                intersects = edges < self%r

            end associate
        end associate
        !&>

    end function circle_t_intersects

end module quad_types
