      subroutine p_gas(rho, u, p, c)
      
c----------------------------------------------------------------------
c   Gamma law EOS: subroutine to calculate the pressure and sound  
 
c     rho    : Density                                              [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
          
      implicit none
      double precision rho, u, p, c   
      double precision gamma 
          
c      For air (idea gas)
c      See Equ.(3.82)

      gamma=1.4
      p = (gamma-1) * rho * u     
      c = sqrt((gamma-1) * u) 
     
      end         
      
      subroutine p_art_water(rho, p, c)
      
c----------------------------------------------------------------------
c   Artificial equation of state for the artificial compressibility 

c     rho    : Density                                              [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
c     Equation of state for artificial compressibility   

      implicit none
      double precision rho, u, p, c
      double precision gamma, rho0

c     Artificial EOS, Form 1 (Monaghan, 1994) 
c     See Equ.(4.88)
c      gamma=7.
c      rho0=1000.       
c      b = 1.013e5
c      p = b*((rho/rho0)**gamma-1)      
c      c = 1480.

c     Artificial EOS, Form 2 (Morris, 1997)
c     See Equ.(4.89)
      c = 0.01
      p = c**2 * rho      
      
      end 