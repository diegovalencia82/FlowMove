
      subroutine viscous_force_mud(mspace,ntype,npairs,pairs,nfilas,dwdx &
          ,epsilon,viscforce,mud)

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

!c     mud(i) : Modelo de viscosidad                      [in]
      
      implicit none
      include 'param.inc'

      integer i,j,npairs,ntype(2)
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1))
      double precision viscforce(9,ntype(1)),epsilon(6,ntype(1))
      double precision sumvisxx,sumvisxy,sumvisxz,sumvisyx,sumvisyy, &
           sumvisyz,sumviszx,sumviszy,sumviszz,eirhoixx,eirhoixy,   &
           eirhoixz,eirhoiyx,eirhoiyy,eirhoiyz,eirhoizx,eirhoizy,   &
           eirhoizz,rhoi2,rhoj2
      double precision mud(ntype(1))

!c      mu = 1.!80.8                   !1.0e-3               ! Viscosidad dinámica
!c      mu = 1.0e-3               ! Viscosidad dinámica

      
      do i=1,ntype(1)-ntype(2)

         sumvisxx = 0.0
         sumvisxy = 0.0
         sumvisxz = 0.0
         sumvisyx = 0.0
         sumvisyy = 0.0
         sumvisyz = 0.0
         sumviszx = 0.0
         sumviszy = 0.0
         sumviszz = 0.0

         rhoi2 = mspace(9,i)**2
         eirhoixx = mud(i) * epsilon(1,i)/rhoi2
         eirhoixy = mud(i) * epsilon(2,i)/rhoi2
         eirhoixz = mud(i) * epsilon(3,i)/rhoi2
         eirhoiyx = eirhoixy
         eirhoiyy = mud(i) * epsilon(4,i)/rhoi2
         eirhoiyz = mud(i) * epsilon(5,i)/rhoi2
         eirhoizx = eirhoixz
         eirhoizy = eirhoiyz
         eirhoizz = mud(i) * epsilon(6,i)/rhoi2
         
         
         do j=1,nfilas(i)
            rhoj2 = mspace(9,pairs(j,i))*mspace(9,pairs(j,i))
            
            sumvisxx = sumvisxx + mspace(8,pairs(j,i)) * ( eirhoixx +  &
            mud(pairs(j,i)) * epsilon(1,pairs(j,i))/rhoj2 ) *dwdx(1,j,i)
            
            sumvisxy = sumvisxy + mspace(8,pairs(j,i)) * ( eirhoixy + &
            mud(pairs(j,i)) * epsilon(2,pairs(j,i))/rhoj2 ) *dwdx(2,j,i)
            
            sumvisxz = sumvisxz + mspace(8,pairs(j,i)) * ( eirhoixz + &
            mud(pairs(j,i)) * epsilon(3,pairs(j,i))/rhoj2 ) *dwdx(3,j,i)

            sumvisyx = sumvisyx + mspace(8,pairs(j,i)) * ( eirhoiyx + &
            mud(pairs(j,i)) * epsilon(2,pairs(j,i))/rhoj2 ) *dwdx(1,j,i)
            
            sumvisyy = sumvisyy + mspace(8,pairs(j,i)) * ( eirhoiyy + &
            mud(pairs(j,i)) * epsilon(4,pairs(j,i))/rhoj2 ) *dwdx(2,j,i)
            
            sumvisyz = sumvisyz + mspace(8,pairs(j,i)) * ( eirhoiyz + &
            mud(pairs(j,i)) * epsilon(5,pairs(j,i))/rhoj2 ) *dwdx(3,j,i)
            
            sumviszx = sumviszx + mspace(8,pairs(j,i)) * ( eirhoizx + &
            mud(pairs(j,i)) * epsilon(3,pairs(j,i))/rhoj2 ) *dwdx(1,j,i)

            sumviszy = sumviszy + mspace(8,pairs(j,i)) * ( eirhoizy + &
            mud(pairs(j,i)) * epsilon(5,pairs(j,i))/rhoj2 ) *dwdx(2,j,i)

            sumviszz = sumviszz + mspace(8,pairs(j,i)) * ( eirhoizz + &
            mud(pairs(j,i)) * epsilon(6,pairs(j,i))/rhoj2 ) *dwdx(3,j,i)            
         enddo

         viscforce(1,i) = muc * sumvisxx
         viscforce(2,i) = muc * sumvisxy
         viscforce(3,i) = muc * sumvisxz
         viscforce(4,i) = muc * sumvisyx
         viscforce(5,i) = muc * sumvisyy
         viscforce(6,i) = muc * sumvisyz
         viscforce(7,i) = muc * sumviszx
         viscforce(8,i) = muc * sumviszy
         viscforce(9,i) = muc * sumviszz
         
      enddo
      
    end subroutine viscous_force_mud
