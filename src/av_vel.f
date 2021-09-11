      subroutine av_vel(ntotal,mass,niac,pair_i,pair_j,
     &           w, vx, rho, av)

c----------------------------------------------------------------------
c     Subroutine to calculate the average velocity to correct velocity
c     for preventing.penetration (monaghan, 1992)

c     ntotal : Number of particles                                  [in]
c     mass   : Particle masses                                      [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     w      : Kernel for all interaction pairs                     [in]
c     vx     : Velocity of each particle                            [in]
c     rho    : Density of each particle                             [in]
c     av     : Average velocityof each particle                    [out]
    
      implicit none
      include 'param.inc'
      
      integer ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction)
      double precision   mass(maxn),w(max_interaction),
     &       vx(dim,maxn), rho(maxn), av(dim, maxn)       
      integer i,j,k,d       
      double precision   vcc, dvx(dim), epsilon
      
c     epsilon --- a small constants chosen by experience, may lead to instability.
c     for example, for the 1 dimensional shock tube problem, the E <= 0.3

      epsilon = 0.3
      
      do i = 1, ntotal
        do d = 1, dim
          av(d,i) = 0.
        enddo 
      enddo
     
      do k=1,niac       
        i = pair_i(k)
        j = pair_j(k)       
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j)            
          av(d, i) = av(d,i) - 2*mass(j)*dvx(d)/(rho(i)+rho(j))*w(k)
          av(d, j) = av(d,j) + 2*mass(i)*dvx(d)/(rho(i)+rho(j))*w(k)                      
        enddo                    
      enddo  
        
      do i = 1, ntotal
        do d = 1, dim
          av(d,i) = epsilon * av(d,i)
        enddo 
      enddo             

      end