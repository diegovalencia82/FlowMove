
      subroutine sph_presion(mspace,ntype,npairs,pairs,nfilas,w,dwdx)

!c     mspace(1,i) = id -- label of each particle                   [out]      
!c     mspace(2 to 4,i) = x-- coordinates of particles              [in/out]
!c     mspace(5 to 7,i) = vx-- velocities of particles              [in/out]
!c     mspace(8,i) = mass-- mass of particles                       [in]
!c     mspace(9,i) = rho-- dnesities of particles                   [in/out]
!c     mspace(10,i) = p-- pressure of particles                     [in/out]
!c     mspace(11,i) = u-- internal energy of particles              [in/out]
!c     mspace(12,i) = itype-- types of particles                    [in]
!c     mspace(13,i) = hsml-- smoothing lengths of particles         [in/out]
!c     mspace(14,i) = c-- sound velocity of particle                [out]
!c     mspace(15,i) = s-- entropy of particles                      [out]
!c     mspace(16,i) = e-- total energy of particles                 [out]
!c     mspace(17 to 19,i) = dx  : dx = vx = dx/dt                   [out]
!c     mspace(20 to 22,i) = dvx = dvx/dt, force per unit mass       [out]
!c     mspace(23,i) = du        : du = du/dt                        [out]
!c     mspace(24,i) = ds        : ds = ds/dt                        [out]
!c     mspace(25,i) = drho      : drho = drh,o/dt                   [out]
!c     mspace(26,i) = eta_c     : Coeficiente de Viscosidad No Lineal [in]      
      
!c     npairs    : maximum number of pairs interaction             [in]
!c     nfilas    : number of interaction for each fluid particle   [in]
!c     mrij      : matrix of ri-rj for all fluid particles for each interaction [in]
!c     mrij      : matrix of xi-xj for all fluid particles for each interaction [in]                
!c     w         : kernel for all interaction pairs                     [in]
!c     dwdx      : Derivative of kernel with respect to x, y and z      [in]

      implicit none
      include 'param.inc'

      integer i,j,npairs,ntype(2)
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      double precision t(fintime+1),dt
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1)), &
          mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))
!c      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))
      double precision sumx,sumy,sumz,pirhoi,pjrhoj,mprhoij
      
! mspace(5,i) = vx
! mspace(6,i) = vy
! mspace(7,i) = vz
! mspace(8,i) = mass
! mspace(9,i) = rho      
! mspace(10,i) = p
      
      do i=1,ntype(1)-ntype(2)
         mspace(20,i) = 0.0
         mspace(21,i) = 0.0
         mspace(22,i) = 0.0
         pirhoi = mspace(10,i)/mspace(9,i)**2
         sumx = 0
         sumy = 0
         sumz = 0
!c         write(*,*)i,'sph presion----',mspace(2,i),mspace(4,i)
         do j=1,nfilas(i)
            pjrhoj = mspace(10,pairs(j,i))/mspace(9,pairs(j,i))**2
            mprhoij = mspace(8,pairs(j,i)) * (pirhoi+pjrhoj)
!c            write(*,*)i,j,mspace(2,pairs(j,i)),mspace(4,pairs(j,i))
!c            write(*,*)'i,j',i,j,pairs(j,i),
!c     +           mspace(9,i),mspace(9,pairs(j,i)),
!c     +           mspace(10,i),mspace(10,pairs(j,i)),(pirhoi+pjrhoj),
!c     +           int(mspace(12,i)),int(mspace(12,pairs(j,i))),
!c     +           dwdx(1,j,i),dwdx(3,j,i)
            sumx = sumx + mprhoij*dwdx(1,j,i)
            sumy = sumy + mprhoij*dwdx(2,j,i)
            sumz = sumz + mprhoij*dwdx(3,j,i)
         enddo
         mspace(20,i) = mspace(20,i) - sumx
         mspace(21,i) = mspace(21,i) - sumy
         mspace(22,i) = mspace(22,i) - sumz
         !c        write(*,*)'sph_presion',i,mspace(20,i),mspace(21,i),mspace(22,i)
      enddo
      
    end subroutine sph_presion
