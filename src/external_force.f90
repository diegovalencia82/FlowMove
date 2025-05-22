
      subroutine external_force(mspace,ntype,npairs,pairs,nfilas,mrij &
           ,mxij,extforce)

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
!c     mxij      : matrix of xi-xj for all fluid particles for each interaction [in]

      implicit none
      include 'param.inc'

      integer i,j,npairs,ntype(2)
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))
      double precision extforce(3,ntype(1)),dd,p1,p2,f,fx,fy,fz,r0

      double precision rijv,xijv,yijv,zijv,r3
      double precision reducfactor

      reducfactor = 1.
      r3 = r0 / 3.
      
      do i=1,ntype(1)-ntype(2)
         extforce(1,i) = 0.0d0
         extforce(2,i) = 0.0d0
         extforce(3,i) = 0.0d0
      enddo
      
      do i=1,ntype(1)-ntype(2)
         extforce(3,i) = -g
      enddo

     
    end subroutine external_force

