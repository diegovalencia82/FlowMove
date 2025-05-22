
      subroutine wijdwij(mspace,pairs,mrij,mxij,npairs,nfilas,ntype,w, &
          dwdx)

!c---------------------------------------------------------------------

!c     Subroutine to get smoothing function pairs
!c     mspace(1,i) = id -- label of each particle                   [in]      
!c     mspace(2 to 4,i) = x-- coordinates of particles              [in]
!c     mspace(5 to 7,i) = vx-- velocities of particles              [in]
!c     npairs    : maximum number of pairs interaction             [in]
!c     nfilas    : number of interaction for each fluid particle   [in]
!c     mrij      : matrix of ri-rj for all fluid particles for each interaction [in]
!c     mxij      : matrix of xi-xj for all fluid particles for each interaction [in]                
!c     w         : kernel for all interaction pairs                     [out]
!c     dwdx      : Derivative of kernel with respect to x, y and z      [out]

      
      implicit none
      include 'param.inc'

      integer i,j
      integer npairs,ntype(2),pairs(npairs,ntype(1)),nfilas(ntype(1))
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      double precision mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1))
      double precision dx(dim),r,dwdx0(dim)
      double precision mhsml
     
      do i=1,ntype(1)-ntype(2)
         do j = 1,nfilas(i)
!c            write(*,*)i,j,nfilas(i),mrij(j,i),mxij(1,j,i),mxij(3,j,i)
!c     +           mxij(2,nfilas(i),i),mxij(3,nfilas(i),i),  
!c     +           mvij(nfilas(i),i),mvxij(1,nfilas(i),i),
!c     +           mvxij(2,nfilas(i),i),mvxij(3,nfilas(i),i) 
            r = mrij(j,i)
            if(dim.eq.2)then
               dx(1) = mxij(1,j,i)
               dx(2) = mxij(3,j,i)
            endif
            if(dim.eq.3)then
               dx(1) = mxij(1,j,i)
               dx(2) = mxij(2,j,i)
               dx(3) = mxij(3,j,i)
            endif
!c            mhsml = ( mspace(13,i)+mspace(13,pairs(j,i)) ) / 2.
!c            call kernel(r,dx,mhsml,w(j,i),dwdx0)
            call kernel(r,dx,mspace(13,i),w(j,i),dwdx0)

!c            write(*,*)j,i,pairs(j,i),r,mspace(13,i),w(j,i)
!hsml = mspace(13,i)
            if(dim.eq.2)then
               dwdx(1,j,i) = dwdx0(1)
               dwdx(2,j,i) = 0.0
               dwdx(3,j,i) = dwdx0(2)
!c               write(*,*)i,j,pairs(j,i),w(j,i)
!c               write(*,*)i,j,pairs(j,i),r,w(j,i),dwdx(1,j,i),dwdx(3,j,i)
!c     +              ,dwdx(2,j,i)
            endif
            if(dim.eq.3)then
               dwdx(1,j,i) = dwdx0(1)
               dwdx(2,j,i) = dwdx0(2)
               dwdx(3,j,i) = dwdx0(3)
            endif
         enddo
      enddo

    end subroutine wijdwij
