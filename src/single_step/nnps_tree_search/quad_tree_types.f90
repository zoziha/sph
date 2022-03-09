!> 四叉树建模
module quad

    use quad_types, only: rectangle_t, point_t, circle_t
    use easy_string_m, only: to_string
    use sph_kinds, only: rk
    implicit none

    !> 四叉树
    type quad_tree_t

        type(rectangle_t) :: boundary
        integer :: capacity                             !! SPH粒子法中，该变量应该设为1
        type(point_t), allocatable :: points(:)
        logical :: divided

        !> 采用可分配内存来存储四叉树内容
        type(quad_tree_t), allocatable :: children(:)

    contains

        procedure :: constructor
        procedure :: insert
        procedure :: subdivide
        procedure :: show

        generic :: query => circle_t_query, rectangle_t_query
        procedure :: circle_t_query, rectangle_t_query
        procedure :: delete_tree

    end type quad_tree_t

contains

    !> 构建四叉树
    subroutine constructor(self, boundary, n)

        class(quad_tree_t), intent(inout) :: self
        type(rectangle_t), intent(in) :: boundary
        integer, intent(in) :: n

        self%boundary = boundary
        self%capacity = n
        allocate (self%points(0))
        self%divided = .false.

        if (boundary%w /= boundary%h) then
            write (*, *) "*<ERROR>* heigh /= width!"
            write (*, *) boundary%w, boundary%h
        end if

    end subroutine constructor

    !> 插入点
    recursive logical function insert(self, point) result(done)

        class(quad_tree_t), intent(inout) :: self
        type(point_t), intent(in) :: point

        integer :: i

        done = .false.

        !> 判断粒子是否处于矩形内部
        if (.not. self%boundary%contains(point)) then
            return
        end if

        if (size(self%points) < self%capacity) then

            !> 添加粒子
            self%points = [self%points, point]
            done = .true.

            !> 当前节点的矩形最大粒子容量已满，指向下属节点
        else

            !> 查询是否有下属节点，无则创建！
            if (.not. self%divided) then

                allocate (self%children(4))

                call self%subdivide()

            end if

            !\TODO: 可能无法运行, Fortran语法？(可以使用递归)
            !> 在下属节点依次尝试插入粒子
            do i = 1, size(self%children)

                if (self%children(i)%insert(point)) then
                    done = .true.
                    exit
                end if

            end do

        end if

        if (.not. done) then
            write (*, *) "*<ERROR>* insert failed!"
            write (*, *) point
            write (*, *) self%boundary
            if (self%divided) then
                do i = 1, size(self%children)
                    write (*, *) self%children(i)%boundary
                end do
            end if
            stop
        end if

    end function insert

    !> 细分
    !> 细分算法会导致某些在矩形边上的粒子无法插入树型表
    subroutine subdivide(self)

        class(quad_tree_t), intent(inout) :: self
        type(rectangle_t) :: nw, ne, sw, se

        associate (x => self%boundary%x, &
                   y => self%boundary%y, &
                   w => self%boundary%w, &
                   h => self%boundary%h)

            !> children: 1-ne(东北), 2-nw(西北), 3-se(东南), 4-sw(西南)

            ! 添加容差 0.0005，使得四叉树存在小的重叠区域，以涵盖所有粒子，消除精度误差；
            ! 存在一些粒子在两个相邻的矩形公共边上，计算机误差使得这些粒子无法插入树型表；
            ! 同时使得，被父矩形包含的粒子，在子矩形中必定能插入。
            ne = rectangle_t(x + 0.25_rk*w, y + 0.25_rk*h, 0.5005_rk*w, 0.5005_rk*h)
            nw = rectangle_t(x - 0.25_rk*w, y + 0.25_rk*h, 0.5005_rk*w, 0.5005_rk*h)
            se = rectangle_t(x + 0.25_rk*w, y - 0.25_rk*h, 0.5005_rk*w, 0.5005_rk*h)
            sw = rectangle_t(x - 0.25_rk*w, y - 0.25_rk*h, 0.5005_rk*w, 0.5005_rk*h)

            call self%children(1)%constructor(ne, self%capacity)
            call self%children(2)%constructor(nw, self%capacity)
            call self%children(3)%constructor(se, self%capacity)
            call self%children(4)%constructor(sw, self%capacity)

            self%divided = .true.

        end associate

    end subroutine subdivide

    !> 打印四叉树
    recursive subroutine show(self, indent)

        class(quad_tree_t), intent(inout) :: self
        integer, intent(in) :: indent

        integer :: i

        do i = 1, size(self%points)
            if (allocated(self%points)) then
                write (*, "(2(2x,f6.2),2x,i0)") self%points(i)%x, self%points(i)%y, self%points(i)%index
            end if
        end do

        if (allocated(self%children)) then

            do i = 1, size(self%children)
                call self%children(i)%show(indent + 1)
            end do

        end if

    end subroutine show

    !> 析构例程（没实际意义）
    subroutine delete_tree(self)

        class(quad_tree_t), intent(out) :: self

    end subroutine delete_tree

    !> 粒子查询操作(矩形)
    recursive subroutine rectangle_t_query(self, range, found)

        class(quad_tree_t), intent(inout) :: self
        type(rectangle_t), intent(in) :: range
        type(point_t), intent(inout), allocatable :: found(:)

        integer :: i

        if (.not. allocated(found)) allocate (found(0))

        if (.not. range%intersects(self%boundary)) then
            return
        end if

        do i = 1, size(self%points)
            if (range%contains(self%points(i))) found = [found, self%points(i)]
        end do

        if (self%divided) then
            do i = 1, size(self%children)
                call self%children(i)%query(range, found)
            end do
        end if

    end subroutine rectangle_t_query

    !> 粒子查询操作(圆形)
    recursive subroutine circle_t_query(self, range, found)

        class(quad_tree_t), intent(inout) :: self
        type(circle_t), intent(in) :: range
        type(point_t), intent(inout), allocatable :: found(:)

        integer :: i

        if (.not. range%intersects(self%boundary)) then
            return
        end if

        if (.not. allocated(found)) allocate (found(0))

        do i = 1, size(self%points)
            if (range%contains(self%points(i))) found = [found, self%points(i)]
        end do

        if (self%divided) then
            do i = 1, size(self%children)
                call self%children(i)%query(range, found)
            end do
        end if

    end subroutine circle_t_query

end module quad
