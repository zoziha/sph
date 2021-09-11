      subroutine init_grid(ntotal,hsml,grid,ngridx,ghsmlx,
     &           maxgridx,mingridx,dgeomx)

c----------------------------------------------------------------------      
c   Subroutine to established a pair linked list by sorting grid cell.
c   It is suitable for a homogeneous particle distribution with the 
c   same smoothing length in an instant. A fixed number of particles
c   lie in each cell. 

c     ntotal   : Number of particles                                [in]
c     hsml     : Smoothing Length                                   [in]
c     grid     : array of grid cells                               [out]
c     ngridx   : Number of sorting grid cells in x, y, z-direction [out]
c     ghsmlx   : Smoothing length measured in cells of the grid    [out]
c     maxgridx : Maximum x-, y- and z-coordinate of grid range     [out]
c     mingridx : Minimum x-, y- and z-coordinate of grid range     [out]
c     dgeomx   : x-, y- and z-expansion of grid range              [out]

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
      integer ntotal, grid(maxngx,maxngy,maxngz), ngridx(dim), 
     &        ghsmlx(dim)
      double precision hsml, maxgridx(dim), mingridx(dim), dgeomx(dim)
      integer i, j, k, d, maxng(dim), ngrid(3)
      double precision nppg

c     Averaged number of particles per grid cell

      parameter( nppg = 3.e0 )

c     Initialize parameters: Maximum number of grid cells

      maxng(1) = maxngx
      if (dim.ge.2) then
        maxng(2) = maxngy
        if (dim.eq.3) then
          maxng(3) = maxngz
        endif
      endif
      
      do d=1,3
        ngrid(d) = 1
      enddo
      
c     Range of sorting grid

      maxgridx(1) = x_maxgeom
      mingridx(1) = x_mingeom
      if (dim.ge.2) then
        maxgridx(2) = y_maxgeom
        mingridx(2) = y_mingeom
        if (dim.eq.3) then
          maxgridx(3) = z_maxgeom
          mingridx(3) = z_mingeom
        endif
      endif

      do d=1,dim
         dgeomx(d) = maxgridx(d) - mingridx(d)
      enddo

c     Number of grid cells in x-, y- and z-direction:

      if (dim.eq.1) then
        ngridx(1) = min(int(ntotal/nppg) + 1,maxng(1))
      else if (dim.eq.2) then
        ngridx(1) = min(
     &  int(sqrt(ntotal*dgeomx(1)/(dgeomx(2)*nppg))) + 1,maxng(1))
        ngridx(2) = min(
     &  int(ngridx(1)*dgeomx(2)/dgeomx(1)) + 1,maxng(2))
      else if (dim.eq.3) then
        ngridx(1) = min(int((ntotal*dgeomx(1)*dgeomx(1)/
     &      (dgeomx(2)*dgeomx(3)*nppg))**(1.e0/3.e0)) + 1,maxng(1))
        ngridx(2) = min(
     &      int(ngridx(1)*dgeomx(2)/dgeomx(1)) + 1,maxng(2))
        ngridx(3) = min(
     &  int(ngridx(1)*dgeomx(3)/dgeomx(1)) + 1,maxng(3))
      endif

c     Smoothing Length measured in grid cells:

      do d=1,dim
         ghsmlx(d) = int(real(ngridx(d))*hsml/dgeomx(d)) + 1
      enddo

      do d=1,dim
        ngrid(d) = ngridx(d)
      enddo

c     Initialize grid

      do i=1,ngrid(1)
        do j=1,ngrid(2)
          do k=1,ngrid(3)
            grid(i,j,k) = 0
          enddo
        enddo
      enddo

      end
