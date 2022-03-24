! 控制台颜色
module console_color_m

    use M_attr, only: attr_ => attr
    use swift_file_m, only: is_windows
    implicit none
    private
    
    public :: attr, is_intel_windows
    
contains

    function attr(str) result(r)
        character(*), intent(in) :: str
        character(:), allocatable :: r
        integer c
        
        if (is_intel_windows()) then
            c = index(str, '>')
            r = str(c+1:)
        else
            r = attr_(str)
        endif
        
    end function attr
    
    function is_intel_windows() result(r)
        logical :: r
#ifdef __INTEL_COMPILER
        r = is_windows()
#else
        r = .false.
#endif
    end function is_intel_windows

end module console_color_m
