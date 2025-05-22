
      subroutine momento(mspace,ntype,npairs,pairs,nfilas,w,dwdx,mrij, &
          mxij,mvij,mvxij,pairsv,nfilasv,mrijv,mxijv)
      
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
!c     mvij      : matrix of vi-vj for all fluid particles for each interaction [in]
!c     mvxij     : matrix of vxi-vxj for all fluid particles for each interaction [in]
      
!c     w         : kernel for all interaction pairs                     [in]
!c     dwdx      : Derivative of kernel with respect to x, y and z      [in]

!c     epsilon(1,i) : Tasa de deformación $\varepsilon^{xx}$
!c     epsilon(2,i) : Tasa de deformación $\varepsilon^{xy} = \varepsilon^{yx}$
!c     epsilon(3,i) : Tasa de deformación $\varepsilon^{xz} = \varepsilon^{zx}$
!c     epsilon(4,i) : Tasa de deformación $\varepsilon^{yy}$
!c     epsilon(5,i) : Tasa de deformación $\varepsilon^{yz} = \varepsilon^{zy}$      
!c     epsilon(6,i) : Tasa de deformación $\varepsilon^{zz}$      
      
      implicit none
      include 'param.inc'

      integer i,j,npairs,ntype(2)
      double precision, dimension(26, nmax) :: mspace
      double precision t(fintime+1),dt
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      integer pairsv(npairs,ntype(1)),nfilasv(ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1)),  &
           mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))            &
           ,mrijv(npairs,ntype(1)),mxijv(3,npairs,ntype(1))
      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))
      double precision epsilon(6,ntype(1)),viscforce(9,ntype(1))
      double precision extforce(3,ntype(1))
      double precision epsilonq(ntype(1)),mud(ntype(1))

      do i=1,ntype(1)!ntype(1)-ntype(2)
         mspace(20,i)=0.
         mspace(21,i)=0.
         mspace(22,i)=0.
      enddo
      
      call sph_presion(mspace,ntype,npairs,pairs,nfilas,w,dwdx) ! 2D y 3D

      epsilon = 0.0

      if ( mvisc == 1)then
         call tasa_deformacion_epsilon(mspace,ntype,npairs,pairs,nfilas,w, &
              dwdx,mxij,mvij,mvxij,epsilon) ! 2D y 3D

         viscforce = 0.0
         call viscous_force(mspace,ntype,npairs,pairs,nfilas,dwdx,epsilon, &
              viscforce)         
      endif

      if ( mvisc == 2)then
         epsilonq = 0.0
         call tasa_deformacion_epsilon2(mspace,ntype,npairs,pairs,nfilas,w, &
              dwdx,mxij,mvij,mvxij,epsilon,epsilonq) ! 2D y 3D
      
!c      do i=1,ntype(1)
!c         write(*,*)i,epsilonq(i),epsilon(1,i),epsilon(2,i),epsilon(3,i),
!c     +        epsilon(4,i),epsilon(5,i),epsilon(6,i)
!c      enddo

         call viscous_model(mspace,ntype,epsilonq,mud) ! 2D y 3D
         
         viscforce = 0.0

         call viscous_force_mud(mspace,ntype,npairs,pairs,nfilas,dwdx, &
              epsilon,viscforce,mud) ! 2D y 3D
      endif
      

      call external_force(mspace,ntype,npairs,pairs,nfilas,mrij,mxij &
          ,extforce)
      
      do i=1,ntype(1) - ntype(2)
!c         write(*,*)'xx',i,mspace(20,i),viscforce(1,i) + viscforce(2,i) +
!c     +        viscforce(3,i)
!c         write(*,*)'yy',i,mspace(21,i),viscforce(4,i) + viscforce(5,i) +
!c     +        viscforce(6,i)
!c         write(*,*)'zz',i,mspace(22,i),viscforce(7,i) + viscforce(8,i) +
!c     +        viscforce(9,i)
         
         mspace(20,i) = mspace(20,i) + viscforce(1,i) + viscforce(2,i) + & 
              viscforce(3,i) + extforce(1,i)
         mspace(21,i) = mspace(21,i) + viscforce(4,i) + viscforce(5,i) + &
              viscforce(6,i) + extforce(2,i)
         mspace(22,i) = mspace(22,i) + viscforce(7,i) + viscforce(8,i) + &
              viscforce(9,i) + extforce(3,i)
      enddo

!c      do i=1,ntype(1)                  !ntype(1)
!c         do j=1,nfilas(i)
!c            write(*,*)'mm',i,mspace(2,pairs(nfilas(j),i))
!c     +           ,mspace(4,pairs(nfilas(j),i))
!c         write(*,*),mspace(20,i),mspace(22,i)
!c         write(*,*),extforce(1,i),extforce(2,i),extforce(3,i)
!c         enddo
!c      enddo
      
    end subroutine momento
      
