module test_eos_m

    use config_m, only: rk, B, rho0, eos_form, c0 => c
    use eos_m, only: p_gas, p_art_water
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_eos

contains

    subroutine collect_eos(test_suite)
        type(unittest_type), allocatable, intent(out) :: test_suite(:)

        test_suite = [ &
                     new_unittest("sub: p_gas", test_p_gas), &
                     new_unittest("sub: p_art_water", test_p_art_water) &
                     ]
    end subroutine collect_eos

    subroutine test_p_gas(error)
        type(error_type), allocatable, intent(out) :: error
        real(rk) :: p, c
        call p_gas(1000.0_rk, 1.0_rk, p, c)
        call check(error, p, 399.99999999999989_rk, "p_gas%p")
        if (allocated(error)) return
        call check(error, c, 0.63245553203367577_rk, "p_gas%c")

    end subroutine test_p_gas

    subroutine test_p_art_water(error)
        type(error_type), allocatable, intent(out) :: error
        real(rk) :: p, c

        B = 142.8_rk; rho0 = 1000.0_rk; eos_form = 1
        call p_art_water(1010.0_rk, p, c)
        call check(error, p, 10.300928280881005_rk, "p_art_water%p1")
        if (allocated(error)) return
        call check(error, c, 1.0300949191898578_rk, "p_art_water%c1")
        if (allocated(error)) return

        eos_form = 2; c0 = 0.01_rk
        call p_art_water(1000.0_rk, p, c)
        call check(error, p, 0.10000000000000001_rk, "p_art_water%p2")
        if (allocated(error)) return
        call check(error, c, 1.0000000000000000E-002_rk, "p_art_water%c2")
        if (allocated(error)) return
        
        eos_form = 3; 
        call p_art_water(1010.0_rk, p, c)
        call check(error, p, 1.0000000000000000E-003_rk, "p_art_water%p3")
        if (allocated(error)) return
        call check(error, c, 1.0000000000000000E-002_rk, "p_art_water%c3")

    end subroutine test_p_art_water

end module test_eos_m
