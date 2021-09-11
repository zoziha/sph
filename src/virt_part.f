      subroutine virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
     &           rho,u,p,itype) 

c----------------------------------------------------------------------
c   Subroutine to determine the information of virtual particles
c   Here only the Monaghan type virtual particles for the 2D shear
c   cavity driven problem are generated.
c     itimestep : Current time step                                 [in]
c     ntotal : Number of particles                                  [in]
c     nvirt  : Number of virtual particles                         [out]
c     hsml   : Smoothing Length                                 [in|out]
c     mass   : Particle masses                                  [in|out]
c     x      : Coordinates of all particles                     [in|out]
c     vx     : Velocities of all particles                      [in|out]
c     rho    : Density                                          [in|out]
c     u      : internal energy                                  [in|out]
c     itype   : type of particles                               [in|out]

      implicit none
      include 'param.inc'
      integer itimestep, ntotal, nvirt, itype(maxn)
      double precision hsml(maxn),mass(maxn),x(dim,maxn),vx(dim,maxn),
     &                 rho(maxn), u(maxn), p(maxn)
      integer i, j, d, im, mp
      double precision xl, dx, v_inf

      if (vp_input) then          
                        
        open(1,file="./test/data/xv_vp.dat")
        open(2,file="./test/data/state_vp.dat")
        open(3,file="./test/data/other_vp.dat")            
        read(1,*) nvirt
        do j = 1, nvirt   
          i = ntotal + j      
          read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
          read(2,*)im, mass(i), rho(i), p(i), u(i)        
          read(3,*)im, itype(i), hsml(i)                            
        enddo  
        close(1)
        close(2) 
        close(3) 
      
	else 
       
	nvirt = 0
        mp = 40
	xl = 1.0e-3
	dx = xl / mp
	v_inf = 1.e-3

c     Monaghan type virtual particle on the Upper side

        do i = 1, 2*mp+1
   	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = (i-1)*dx/2 
          x(2, ntotal + nvirt) = xl  
          vx(1, ntotal + nvirt) = v_inf
	  vx(2, ntotal + nvirt) = 0.
        enddo

c     Monaghan type virtual particle on the Lower side

        do i = 1, 2*mp+1
   	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = (i-1)*dx/2 
          x(2, ntotal + nvirt) = 0.  
          vx(1, ntotal + nvirt) = 0.
	  vx(2, ntotal + nvirt) = 0.
        enddo

c     Monaghan type virtual particle on the Left side

        do i = 1, 2*mp-1
   	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = 0. 
          x(2, ntotal + nvirt) = i*dx/2
          vx(1, ntotal + nvirt) = 0.
	  vx(2, ntotal + nvirt) = 0.
        enddo

c     Monaghan type virtual particle on the Right side

        do i = 1, 2*mp-1
   	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = xl 
          x(2, ntotal + nvirt) = i*dx/2  
          vx(1, ntotal + nvirt) = 0.
	  vx(2, ntotal + nvirt) = 0.
        enddo

	do i = 1, nvirt
	  rho (ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
	  p(ntotal + i) = 0.
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -2
	  hsml(ntotal + i) = dx
        enddo
        
      endif   

      if (mod(itimestep,save_step).eq.0) then
        open(1,file="./test/data/xv_vp.dat")
        open(2,file="./test/data/state_vp.dat")
        open(3,file="./test/data/other_vp.dat")            
        write(1,*) nvirt
        do i = ntotal + 1, ntotal + nvirt         
          write(1,1001) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)              
          write(2,1002) i, mass(i), rho(i), p(i), u(i)
          write(3,1003) i, itype(i), hsml(i)                               
        enddo       
1001    format(1x, I6, 6(2x, e14.8))
1002    format(1x, I6, 7(2x, e14.8)) 
1003    format(1x, I6, 2x, I4, 2x, e14.8)
        close(1)
        close(2) 
        close(3) 
      endif 

      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
         print *,' >> Statistics: Virtual boundary particles:'
         print *,'          Number of virtual particles:',NVIRT
        endif     
      endif

      end