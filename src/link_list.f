      subroutine link_list(itimestep, ntotal,hsml,x,niac,pair_i,
     &           pair_j,w,dwdx,countiac)
     
c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing funciton for each particle and
c   the interaction parameters used by the SPH algorithm. Interaction 
c   pairs are determined by using a sorting grid linked list  

c     itimestep : Current time step                                 [in]
c     ntotal    : Number of particles                               [in]
c     hsml      : Smoothing Length, same for all particles          [in]
c     x         : Coordinates of all particles                      [in]
c     niac      : Number of interaction pairs                      [out]
c     pair_i    : List of first partner of interaction pair        [out]
c     pair_j    : List of second partner of interaction pair       [out]
c     w         : Kernel for all interaction pairs                 [out]
c     dwdx      : Derivative of kernel with respect to x, y and z  [out]
c     countiac  : Number of neighboring particles                  [out]

      implicit none
      include 'param.inc'

c     Parameter used for sorting grid cells in the link list algorithm
c     maxngx  : Maximum number of sorting grid cells in x-direction
c     maxngy  : Maximum number of sorting grid cells in y-direction
c     maxngz  : Maximum number of sorting grid cells in z-direction
c     Determining maximum number of sorting grid cells:
c     (For an homogeneous particle distribution:)
c     1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
c     2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
c     3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)
      integer maxngx,maxngy,maxngz
      parameter ( maxngx  = 100        ,
     &            maxngy  = 100        ,
     &            maxngz  = 1          )      
      integer itimestep, ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), countiac(maxn)
      double precision hsml, x(dim,maxn),w(max_interaction),
     &       dwdx(dim,max_interaction)
      integer i, j, d, scale_k, sumiac, maxiac, noiac, miniac, maxp,minp    
      integer grid(maxngx,maxngy,maxngz),xgcell(3,maxn),gcell(3),
     &     xcell,ycell,zcell,celldata(maxn),minxcell(3),maxxcell(3),
     &     dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim)
      double precision hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),
     &       tdwdx(dim), dgeomx(dim)

      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
         scale_k = 3 
      endif 
     
      do i=1,ntotal
        countiac(i) = 0
      enddo

c     Initialize grid:  

      call init_grid(ntotal,hsml,grid,ngridx,ghsmlx,
     &     maxgridx,mingridx,dgeomx)
      
c     Position particles on grid and create linked list:
      
      do i=1,ntotal
        call grid_geom(i,x(1,i),ngridx,maxgridx,mingridx,dgeomx,gcell)
        do d=1,dim
          xgcell(d,i) = gcell(d)
        enddo
        celldata(i) = grid(gcell(1),gcell(2),gcell(3))
        grid(gcell(1),gcell(2),gcell(3)) = i
      enddo

c     Determine interaction parameters:

      niac = 0
      do i=1,ntotal-1

c     Determine range of grid to go through:
         
        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,dim
          dnxgcell(d) = xgcell(d,i) - ghsmlx(d)
          dpxgcell(d) = xgcell(d,i) + ghsmlx(d)
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),ngridx(d))
        enddo

c     Search grid:
      
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = grid(xcell,ycell,zcell)
 1            if (j.gt.i) then
                dx(1) = x(1,i) - x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,dim
                  dx(d) = x(d,i) - x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                  if (niac.lt.max_interaction) then

c     Neighboring pair list, and totalinteraction number and
c     the interaction number for each particle 

                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    r = sqrt(dr)
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1
                           
C--- Kernel and derivations of kernel

                    call kernel(r,dx,hsml,w(niac),tdwdx)
	            do d = 1, dim
	              dwdx(d,niac)=tdwdx(d)
                    enddo                  
                  else
                    print *,
     &              ' >>> Error <<< : too many interactions'
                    stop
                  endif
                endif
                j = celldata(j)
                goto 1
              endif
            enddo
          enddo
        enddo
      enddo

c     Statistics for the interaction

      sumiac = 0
      maxiac = 0
      miniac = 1000
      noiac  = 0
      do i=1,ntotal
        sumiac = sumiac + countiac(i)
        if (countiac(i).gt.maxiac) then
	  maxiac = countiac(i)
	  maxp = i
	endif
	if (countiac(i).lt.miniac) then 
	  miniac = countiac(i)
          minp = i
	endif
        if (countiac(i).eq.0)      noiac  = noiac + 1
      enddo
 
      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
          print *,' >> Statistics: interactions per particle:'
          print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
          print *,'**** Particle:',minp, ' minimal interactions:',miniac
          print *,'**** Average :',real(sumiac)/real(ntotal)
          print *,'**** Total pairs : ',niac
          print *,'**** Particles with no interactions:',noiac
        endif     
      endif

      end
