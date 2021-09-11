      subroutine direct_find(itimestep, ntotal,hsml,x,niac,pair_i,
     &           pair_j,w,dwdx,countiac)

c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing funciton for each particle and
c   the interaction parameters used by the SPH algorithm. Interaction 
c   pairs are determined by directly comparing the particle distance 
c   with the corresponding smoothing length.
c   See p.148 in Chapter 4

c     itimestep : Current time step                                 [in]
c     ntotal    : Number of particles                               [in]
c     hsml      : Smoothing Length                                  [in]
c     x         : Coordinates of all particles                      [in]
c     niac      : Number of interaction pairs                      [out]
c     pair_i    : List of first partner of interaction pair        [out]
c     pair_j    : List of second partner of interaction pair       [out]
c     w         : Kernel for all interaction pairs                 [out]
c     dwdx      : Derivative of kernel with respect to x, y and z  [out]
c     countiac  : Number of neighboring particles                  [out]

      implicit none
      include 'param.inc'
      
      integer itimestep, ntotal,niac,pair_i(max_interaction),
     &        pair_j(max_interaction), countiac(maxn)
      double precision hsml(maxn), x(dim,maxn), w(max_interaction),
     &       dwdx(dim,max_interaction)
      integer i, j, d,  sumiac, maxiac, miniac, noiac,
     &        maxp, minp, scale_k 
      double precision dxiac(dim), driac, r, mhsml, tdwdx(dim)     
c     Smoothing kernel function 
c     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c         = 2, Gauss kernel   (Gingold and Monaghan 1981) 
c         = 3, Quintic kernel (Morris 1997)
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
      
      niac = 0

      do i=1,ntotal-1     
        do j = i+1, ntotal
          dxiac(1) = x(1,i) - x(1,j)
          driac    = dxiac(1)*dxiac(1)
          do d=2,dim
            dxiac(d) = x(d,i) - x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (hsml(i)+hsml(j))/2.
          if (sqrt(driac).lt.scale_k*mhsml) then
            if (niac.lt.max_interaction) then    

c     Neighboring pair list, and totalinteraction number and
c     the interaction number for each particle 

              niac = niac + 1
              pair_i(niac) = i
              pair_j(niac) = j
              r = sqrt(driac)
              countiac(i) = countiac(i) + 1
              countiac(j) = countiac(j) + 1

c     Kernel and derivations of kernel

              call kernel(r,dxiac,mhsml,w(niac),tdwdx)
              do d=1,dim
                dwdx(d,niac) = tdwdx(d)
              enddo                                  	     
            else
              print *,
     &        ' >>> ERROR <<< : Too many interactions' 
              stop
            endif
          endif
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