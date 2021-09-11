      subroutine grid_geom(i,x,ngridx,maxgridx,mingridx,dgeomx,xgcell)

c----------------------------------------------------------------------
c   Subroutine to calculate the coordinates (xgcell) of the cell of 
c   the sorting  grid, in which the particle with coordinates (x) lies.

c     x        : Coordinates of particle                            [in]    
c     ngridx   : Number of sorting grid cells in x, y, z-direction  [in]
c     maxgridx : Maximum x-, y- and z-coordinate of grid range      [in]
c     mingridx : Minimum x-, y- and z-coordinate of grid range      [in]
c     dgeomx   : x-, y- and z-expansion of grid range               [in]
c     xgcell   : x-, y- and z-coordinte of sorting grid cell       [out]

      implicit none
      include 'param.inc'

      integer i, ngridx(dim),xgcell(3)
      double precision x(dim), maxgridx(dim), mingridx(dim), dgeomx(dim)
      integer d

      do d=1,3
        xgcell(d) = 1
      enddo

      do d=1,dim
        if ((x(d).gt.maxgridx(d)).or.(x(d).lt.mingridx(d))) then
          print *,' >>> ERROR <<< : Particle out of range'
          print *,'    Particle position: x(',i,d,') = ',x(d)
          print *,'    Range: [xmin,xmax](',D,') = 
     &         [',mingridx(d),',',maxgridx(d),']'
          stop
        else
          xgcell(d) = int(real(ngridx(d))/dgeomx(d)*
     &         (x(d)-mingridx(d)) + 1.e0)
        endif
      enddo

      end