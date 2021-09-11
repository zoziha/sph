      subroutine time_integration(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt )
     
c----------------------------------------------------------------------
c      x-- coordinates of particles                       [input/output]
c      vx-- velocities of particles                       [input/output]
c      mass-- mass of particles                                  [input]
c      rho-- dnesities of particles                       [input/output]
c      p-- pressure  of particles                         [input/output]
c      u-- internal energy of particles                   [input/output]
c      c-- sound velocity of particles                          [output]
c      s-- entropy of particles, not used here                  [output]
c      e-- total energy of particles                            [output]
c      itype-- types of particles                               [input]
c           =1   ideal gas
c           =2   water
c           =3   tnt
c      hsml-- smoothing lengths of particles              [input/output]
c      ntotal-- total particle number                            [input]  
c      maxtimestep-- maximum timesteps                           [input]
c      dt-- timestep                                             [input]
   

      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal, maxtimestep
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn), 
     &       rho(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), 
     &       hsml(maxn), dt
      integer i, j, k, itimestep, d, current_ts, nstart        
      double precision  x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),
     &       rho_min(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn),  
     &       drho(maxn),  av(dim, maxn), ds(maxn), 
     &       t(maxn), tdsdt(maxn)         
      double precision  time, temp_rho, temp_u
               
      do i = 1, ntotal
        do d = 1, dim
          av(d, i) = 0.
        enddo
      enddo  
     
      do itimestep = nstart+1, nstart+maxtimestep   
	   
        current_ts=current_ts+1
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif      
      
c     If not first time step, then update thermal energy, density and 
c     velocity half a time step  

        if (itimestep .ne. 1) then

          do i = 1, ntotal        
            u_min(i) = u(i)
            temp_u=0.
	    if (dim.eq.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)
            u(i) = u(i) + (dt/2.)* (du(i)+temp_u)
            if(u(i).lt.0)  u(i) = 0.                 
            
            if (.not.summation_density) then    
              rho_min(i) = rho(i)
	      temp_rho=0.
	      if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
              rho(i) = rho(i) +(dt/2.)*( drho(i)+ temp_rho)
            endif 
           
            do d = 1, dim
              v_min(d, i) = vx(d, i)
              vx(d, i) = vx(d, i) + (dt/2.)*dvx(d, i)
            enddo
          enddo 
          
        endif

c---  Definition of variables out of the function vector:    
      
        call single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u, s, 
     &       rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av)  
                  
        if (itimestep .eq. 1) then
       
          do i=1,ntotal
            temp_u=0.
   	    if (dim.eq.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)        
            u(i) = u(i) + (dt/2.)*(du(i) + temp_u)
            if(u(i).lt.0)  u(i) = 0.             
         
            if (.not.summation_density ) then
	      temp_rho=0.
	      if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
              rho(i) = rho(i) + (dt/2.)* (drho(i)+temp_rho)
            endif
         
            do d = 1, dim        
              vx(d, i) = vx(d, i) + (dt/2.) * dvx(d, i) + av(d, i)
              x(d, i) = x(d, i) + dt * vx(d, i)		  
            enddo  	            
          enddo 
                  
        else   
                    
          do i=1,ntotal            
            temp_u=0.
	    if (dim.eq.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)                       
            u(i) = u_min(i) + dt*(du(i)+temp_u)
            if(u(i).lt.0)  u(i) = 0.          
            
            if (.not.summation_density ) then 
              temp_rho=0.
	      if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)        	           
              rho(i) = rho_min(i) + dt*(drho(i)+temp_rho)
            endif
                  
            do d = 1, dim                   
              vx(d, i) = v_min(d, i) + dt * dvx(d, i) + av(d, i)
              x(d, i) = x(d, i) + dt * vx(d, i)                  
            enddo
          enddo
        
        endif 
     
        time = time + dt

	if (mod(itimestep,save_step).eq.0) then
          call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
	endif 

        if (mod(itimestep,print_step).eq.0) then
          write(*,*)
          write(*,101)'x','velocity', 'dvx'    
          write(*,100)x(1,moni_particle), vx(1,moni_particle), 
     &                dvx(1,moni_particle)    
        endif
        
101     format(1x,3(2x,a12))	 
100     format(1x,3(2x,e12.6))
	 
      enddo

      nstart=current_ts

      end
