      subroutine sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,
     &           itype,rho)

C----------------------------------------------------------------------
C   Subroutine to calculate the density with SPH summation algorithm.
c   See Equ.(4.35)

C     ntotal : Number of particles                                  [in]
C     hsml   : Smoothing Length                                     [in]
C     mass   : Particle masses                                      [in]
C     niac   : Number of interaction pairs                          [in]
C     pair_i : List of first partner of interaction pair            [in]
C     pair_j : List of second partner of interaction pair           [in]
C     w      : Kernel for all interaction pairs                     [in]
c     itype   : type of particles                                   [in]
c     x       : Coordinates of all particles                        [in]
c     rho    : Density                                             [out]
    
      implicit none
      include 'param.inc'
      
      integer ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn)  
      double precision hsml(maxn),mass(maxn), w(max_interaction),
     &       rho(maxn) 
      integer i, j, k, d      
      double precision selfdens, hv(dim), r, wi(maxn)     

c     wi(maxn)---integration of the kernel itself
        
      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      
c     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*mass(i)/rho(i)
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        wi(i) = wi(i) + mass(j)/rho(j)*w(k)
        wi(j) = wi(j) + mass(i)/rho(i)*w(k)
      enddo

c     Secondly calculate the rho integration over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        rho(i) = selfdens*mass(i)
      enddo

c     Calculate SPH sum for rho:
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        rho(i) = rho(i) + mass(j)*w(k)
        rho(j) = rho(j) + mass(i)*w(k)
      enddo

c     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
     
      if (nor_density) then 
        do i=1, ntotal
          rho(i)=rho(i)/wi(i)
        enddo
      endif 
      
      end
      
      subroutine con_density(ntotal,mass,niac,pair_i,pair_j,
     &           dwdx,vx,itype,x,rho, drhodt)

c----------------------------------------------------------------------
c     Subroutine to calculate the density with SPH continuity approach.
c     See Equ.(4.34)

c     ntotal : Number of particles                                  [in]
c     mass   : Particle masses                                      [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : derivation of Kernel for all interaction pairs       [in]
c     vx     : Velocities of all particles                          [in]
c     itype   : type of particles                                   [in]
c     x      : Coordinates of all particles                         [in]
c     rho    : Density                                              [in]
c     drhodt : Density change rate of each particle                [out]   

      implicit none
      include 'param.inc'
      
      integer ntotal,niac,pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn)    
      double precision mass(maxn), dwdx(dim, max_interaction),
     &       vx(dim,maxn), x(dim,maxn), rho(maxn), drhodt(maxn)
      integer i,j,k,d    
      double precision    vcc, dvx(dim) 
      
      do i = 1, ntotal
        drhodt(i) = 0.
      enddo
     
      do k=1,niac      
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j) 
        enddo        
        vcc = dvx(1)*dwdx(1,k)        
        do d=2,dim
          vcc = vcc + dvx(d)*dwdx(d,k)
        enddo    
        drhodt(i) = drhodt(i) + mass(j)*vcc
        drhodt(j) = drhodt(j) + mass(i)*vcc       
      enddo    
	 
      end