      subroutine input(x, vx, mass, rho, p, u, itype, hsml, ntotal)
      
c----------------------------------------------------------------------
c     Subroutine for loading or generating initial particle information

c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c     hsml-- smoothing lengths of particles                        [out]
c     ntotal-- total particle number                               [out]

      implicit none     
      include 'param.inc'

      integer itype(maxn), ntotal       
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn), 
     &                 p(maxn), u(maxn), hsml(maxn), rho(maxn)
      integer i, d, im       

c     load initial particle information from external disk file

      if(config_input) then    
                        
        open(1,file="./TEST/data/f_xv.dat")
        open(2,file="./TEST/data/f_state.dat")
        open(3,file="./TEST/data/f_other.dat")        
      
        write(*,*)'  **************************************************'
        write(*,*)'      Loading initial particle configuration...   '       
	  read (1,*) ntotal 
        write(*,*)'      Total number of particles   ', ntotal    	
        write(*,*)'  **************************************************'	
        do i = 1, ntotal         
          read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
          read(2,*)im, mass(i), rho(i), p(i), u(i)        
          read(3,*)im, itype(i), hsml(i)                                        
        enddo          
          
      else 
          
        open(1,file="./TEST/data/ini_xv.dat")
        open(2,file="./TEST/data/ini_state.dat")
        open(3,file="./TEST/data/ini_other.dat") 
       
      if (shocktube) call shock_tube(x, vx, mass, rho, p, u, 
     &                    itype, hsml, ntotal)                

      if (shearcavity) call shear_cavity(x, vx, mass, rho, p, u, 
     &                      itype, hsml, ntotal)
        do i = 1, ntotal 
          write(1,1001) i, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim) 
          write(2,1002) i, mass(i), rho(i), p(i), u(i)         
          write(3,1003) i, itype(i), hsml(i)    
        enddo   
1001    format(1x, I5, 6(2x, e14.8)) 
1002    format(1x, I5, 7(2x, e14.8)) 
1003    format(1x, I5, 2x, I2, 2x, e14.8) 
        write(*,*)'  **************************************************'
        write(*,*)'      Initial particle configuration generated   '       
        write(*,*)'      Total number of particles   ', ntotal    	
        write(*,*)'  **************************************************' 

      endif

      close(1)
      close(2) 
      close(3) 

      end              
       
       
      subroutine shock_tube(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     1 d noh shock tube problem
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =1   ideal gas
c     hsml-- smoothing lengths of particles                        [out]
c     ntotal-- total particle number                               [out]

      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, d
      double precision space_x     

      ntotal=400
      space_x=0.6/80.      
      
      do i=1,ntotal
        mass(i)=0.75/400.
        hsml(i)=0.015
        itype(i)=1
        do d = 1, dim
          x(d,i) = 0. 
          vx(d,i) = 0.
        enddo        
      enddo                
                
      do i=1,320
        x(1,i)=-0.6+space_x/4.*(i-1)
      enddo
      
      do i=320+1,ntotal
        x(1,i)=0.+space_x*(i-320)
      enddo               
                           
      do i=1,ntotal
        if (x(1,i).le.1.e-8) then
          u(i)=2.5
          rho(i)=1.
          p(i)=1.
        endif 
        if (x(1,i).gt.1.e-8)  then
          u(i)=1.795
          rho(i)=0.25
          p(i)=0.1795
        endif        
      enddo        
	               
      end
      
      subroutine shear_cavity(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2 d shear driven cavity probem with Re = 1
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =2   water
c     h-- smoothing lengths of particles                           [out]
c     ntotal-- total particle number                               [out]

      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, j, d, m, n, mp, np, k
      double precision xl, yl, dx, dy

c     Giving mass and smoothing length as well as other data.

      m = 41
      n = 41
      mp = m-1
      np = n-1
      ntotal = mp * np
      xl = 1.e-3
      yl = 1.e-3
      dx = xl/mp
      dy = yl/np

      do i = 1, mp
	do j = 1, np
	  k = j + (i-1)*np
	  x(1, k) = (i-1)*dx + dx/2.
	  x(2, k) = (j-1)*dy + dy/2.
        enddo
      enddo

      do i = 1, mp*np
	vx(1, i) = 0.
	vx(2, i) = 0.      
        rho (i) = 1000.   
        mass(i) = dx*dy*rho(i)  
        p(i)= 0.   
        u(i)=357.1
        itype(i) = 2
        hsml(i) = dx
      enddo  

      end	 
