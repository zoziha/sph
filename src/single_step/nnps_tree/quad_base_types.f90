module quad_types

    use sph_kinds, only: rk
    implicit none

    !> 二维点
    type point_t
        real(rk) :: x  !! 点的横坐标
        real(rk) :: y  !! 点的纵坐标
        integer index !! 粒子索引
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
    logical function rectangle_t_contains(self, point) result(contains)

        class(rectangle_t), intent(in) :: self
        type(point_t), intent(in) :: point

        associate (x => self%x, &
                   y => self%y, &
                   w => self%w, &
                   h => self%h)

            contains = (point%x >= x - 0.5*w) .and. &
                       (point%x <= x + 0.5*w) .and. &
                       (point%y >= y - 0.5*h) .and. &
                       (point%y <= y + 0.5*h)

        end associate

    end function rectangle_t_contains

    !> 查询几何形状是否有交集
    logical function rectangle_t_intersects(self, range) result(intersects)

        class(rectangle_t), intent(in) :: self
        type(rectangle_t), intent(in) :: range

        associate (left => self%x - 0.25*self%w, &
                   right => self%x + 0.25*self%w, &
                   bottom => self%y - 0.25*self%h, &
                   top => self%y + 0.25*self%h, &
                   left_ => range%x - 0.25*range%w, &
                   right_ => range%x + 0.25*range%w, &
                   bottom_ => range%y - 0.25*range%h, &
                   top_ => range%y + 0.25*range%h)

            intersects = (left < right_) .or. &
                         (right > left_) .or. &
                         (bottom < top_) .or. &
                         (top > bottom_)

        end associate

    end function rectangle_t_intersects

    !> 查询粒子是否在圆形内
    logical function circle_t_contains(self, point) result(contains)

        class(circle_t), intent(in) :: self
        type(point_t), intent(in) :: point

        contains = hypot(self%x - point%x, self%y - point%y) <= self%r

    end function circle_t_contains

    !> 查询几何形状是否有交集
    logical function circle_t_intersects(self, range) result(intersects)

        class(circle_t), intent(in) :: self
        type(rectangle_t), intent(in) :: range

        real(rk) :: x_dist, y_dist

        x_dist = abs(range%x - self%x)
        y_dist = abs(range%y - self%y)

        associate (r => self%r, &
                   w => range%w, &
                   h => range%h)
            associate (edges => hypot(x_dist - w, y_dist - h))

                !> no intersection
                if (x_dist > r + w .or. y_dist > r + h) then
                    intersects = .false.
                    return
                end if

                !> intersection within the circle
                if (x_dist <= w .or. y_dist <= h) then
                    intersects = .true.
                    return
                end if

                !> intersection on the edge of the circle
                intersects = edges <= self%r

            end associate
        end associate

    end function circle_t_intersects

end module quad_types
