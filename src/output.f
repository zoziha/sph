      subroutine output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal) 
      
c----------------------------------------------------------------------           
c     Subroutine for saving particle information to external disk file

c     x-- coordinates of particles                                  [in]
c     vx-- velocities of particles                                  [in]
c     mass-- mass of particles                                      [in]
c     rho-- dnesities of particles                                  [in]
c     p-- pressure  of particles                                    [in]
c     u-- internal energy of particles                              [in]
c     c-- sound velocity of particles                               [in]
c     itype-- types of particles                                    [in]
c     hsml-- smoothing lengths of particles                         [in]
c     ntotal-- total particle number                                [in]

      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn), 
     &       rho(maxn),p(maxn), u(maxn), c(maxn), hsml(maxn)
      integer i, d, npart     
      
      open(1,file="./TEST/data/f_xv.dat")
      open(2,file="./TEST/data/f_state.dat")
      open(3,file="./TEST/data/f_other.dat") 
     
      write(1,*) ntotal
      do i = 1, ntotal         
        write(1,1001) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)              
        write(2,1002) i, mass(i), rho(i), p(i), u(i)
        write(3,1003) i, itype(i), hsml(i)                               
      enddo 
      
1001  format(1x, I6, 6(2x, e14.8))
1002  format(1x, I6, 7(2x, e14.8)) 
1003  format(1x, I6, 2x, I4, 2x, e14.8)      
                                        
      close(1)
      close(2) 
      close(3) 
      
      end           