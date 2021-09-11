      subroutine art_visc(ntotal,hsml,mass,x,vx,niac,rho,c,
     &           pair_i,pair_j,w,dwdx,dvxdt,dedt)

c----------------------------------------------------------------------
c     Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
c     See Equ.(4.66) Equ.(4.62)

c     ntotal : Number of particles (including virtual particles)    [in]
c     hsml   : Smoothing Length                                     [in]
c     mass   : Particle masses                                      [in]
c     x      : Coordinates of all particles                         [in]
c     vx     : Velocities of all particles                          [in]
c     niac   : Number of interaction pairs                          [in]
c     rho    : Density                                              [in]
c     c      : Temperature                                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     w      : Kernel for all interaction pairs                     [in]
c     dwdx   : Derivative of kernel with respect to x, y and z      [in]
c     dvxdt  : Acceleration with respect to x, y and z             [out] 
c     dedt   : Change of specific internal energy                  [out]
 
      implicit none
      include 'param.inc'
      
      integer ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction)            
      double precision hsml(maxn), mass(maxn), x(dim,maxn),vx(dim,maxn),
     &       rho(maxn), c(maxn), w(max_interaction),
     &       dwdx(dim,max_interaction), dvxdt(dim,maxn), dedt(maxn)
      integer i,j,k,d
      double precision dx, dvx(dim), alpha, beta, etq, piv,
     &       muv, vr, rr, h, mc, mrho, mhsml
     
c     Parameter for the artificial viscosity:
c     Shear viscosity
      parameter( alpha = 1.e0   )
     
c     Bulk viscosity
      parameter( beta  = 1.e0  ) 
      
c     Parameter to avoid singularities
      parameter( etq   = 0.1e0 )
           
      do i=1,ntotal
        do d=1,dim
          dvxdt(d,i) = 0.e0
        enddo
        dedt(i) = 0.e0
      enddo   
     
c     Calculate SPH sum for artificial viscosity
      
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        mhsml= (hsml(i)+hsml(j))/2.
        vr = 0.e0
        rr = 0.e0
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j)
          dx     =  x(d,i) -  x(d,j)
          vr     = vr + dvx(d)*dx
          rr     = rr + dx*dx
        enddo

c     Artificial viscous force only if v_ij * r_ij < 0

        if (vr.lt.0.e0) then

c     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )
            
          muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)
          
c     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

          mc   = 0.5e0*(c(i) + c(j))
          mrho = 0.5e0*(rho(i) + rho(j))
          piv  = (beta*muv - alpha*mc)*muv/mrho              

c     Calculate SPH sum for artificial viscous force

          do d=1,dim
            h = -piv*dwdx(d,k)
            dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
            dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
            dedt(i) = dedt(i) - mass(j)*dvx(d)*h
            dedt(j) = dedt(j) - mass(i)*dvx(d)*h
          enddo
        endif
      enddo

c     Change of specific internal energy:

      do i=1,ntotal
         dedt(i) = 0.5e0*dedt(i)        
      enddo

      end