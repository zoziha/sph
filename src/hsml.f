      subroutine h_upgrade(dt,ntotal, mass, vx, rho, niac, 
     &           pair_i, pair_j, dwdx, hsml)

c-----------------------------------------------------------------------
c     Subroutine to evolve smoothing length

c     dt     : time step                                            [in]
c     ntotal : Number of particles                                  [in]
c     mass   : Particle masses                                      [in]
c     vx     : Velocities of all particles                          [in]
c     rho    : Density                                              [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : Derivative of kernel with respect to x, y and z      [in]
c     hsml   : Smoothing Length                                 [in/out]
   
      implicit none
      include 'param.inc'
      
      integer ntotal, niac, pair_i(max_interaction), 
     &        pair_j(max_interaction)
      double precision mass(maxn), vx(dim, maxn), rho(maxn),
     &       dwdx(dim, max_interaction), hsml(maxn)     
      integer i,j,k,d
      double precision dt, fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn)     

      if (sle.eq.0 ) then     

c---  Keep smoothing length unchanged. 
     
        return
      
      else if (sle.eq.2) then
      
c---  dh/dt = (-1/dim)*(h/rho)*(drho/dt).

        do i=1,ntotal
          vcc(i) = 0.e0
        enddo
      
        do k=1,niac
          i = pair_i(k)
          j = pair_j(k)
          do d=1,dim
            dvx(d) = vx(d,j) - vx(d,i) 
          enddo
          hvcc = dvx(1)*dwdx(1,k)
          do d=2,dim
            hvcc = hvcc + dvx(d)*dwdx(d,k)
          enddo    
          vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
          vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)         
        enddo  
        
        do i = 1, ntotal
          dhsml(i) = (hsml(i)/dim)*vcc(i)
          hsml(i) = hsml(i) + dt*dhsml(i)      
          if (hsml(i).le.0) hsml(i) = hsml(i) - dt*dhsml(i) 
        enddo
    
      else if(sle.eq.1) then
            
        fac = 2.0
        do i = 1, ntotal          
          hsml(i) = fac * (mass(i)/rho(i))**(1./dim)
        enddo
       
      endif 
       
      end 