
      subroutine viscous_force(mspace,ntype,npairs,pairs,nfilas,dwdx, &
          epsilon,viscforce)

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
      
      implicit none
      include 'param.inc'

      integer i,j,npairs,ntype(2)
      double precision, dimension(26, nmax) :: mspace
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1))
      double precision viscforce(9,ntype(1)),epsilon(6,ntype(1))
      double precision sumvisxx,sumvisxy,sumvisxz,sumvisyx,sumvisyy, &
           sumvisyz,sumviszx,sumviszy,sumviszz,eirhoixx,eirhoixy,    &
           eirhoixz,eirhoiyx,eirhoiyy,eirhoiyz,eirhoizx,eirhoizy,    &
           eirhoizz,rhoi2,rhoj2     

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
         eirhoixx = epsilon(1,i)/rhoi2
         eirhoixy = epsilon(2,i)/rhoi2
         eirhoixz = epsilon(3,i)/rhoi2
         eirhoiyx = eirhoixy
         eirhoiyy = epsilon(4,i)/rhoi2
         eirhoiyz = epsilon(5,i)/rhoi2
         eirhoizx = eirhoixz
         eirhoizy = eirhoiyz
         eirhoizz = epsilon(6,i)/rhoi2
         
         
         do j=1,nfilas(i)
            rhoj2 = mspace(9,pairs(j,i))*mspace(9,pairs(j,i))
            
            sumvisxx = sumvisxx + mspace(8,pairs(j,i)) * ( eirhoixx + &
                 epsilon(1,pairs(j,i))/rhoj2 ) * dwdx(1,j,i)
            
            sumvisxy = sumvisxy + mspace(8,pairs(j,i)) * ( eirhoixy + &
                 epsilon(2,pairs(j,i))/rhoj2 ) * dwdx(2,j,i)
            
            sumvisxz = sumvisxz + mspace(8,pairs(j,i)) * ( eirhoixz + &
                 epsilon(3,pairs(j,i))/rhoj2 ) * dwdx(3,j,i)
            
            sumvisyx = sumvisyx + mspace(8,pairs(j,i)) * ( eirhoiyx + &
                 epsilon(2,pairs(j,i))/rhoj2 ) * dwdx(1,j,i)
            
            sumvisyy = sumvisyy + mspace(8,pairs(j,i)) * ( eirhoiyy + &
                 epsilon(4,pairs(j,i))/rhoj2 ) * dwdx(2,j,i)
            
            sumvisyz = sumvisyz + mspace(8,pairs(j,i)) * ( eirhoiyz + &
                 epsilon(5,pairs(j,i))/rhoj2 ) * dwdx(3,j,i)
            
            sumviszx = sumviszx + mspace(8,pairs(j,i)) * ( eirhoizx + &
                 epsilon(3,pairs(j,i))/rhoj2 ) * dwdx(1,j,i)
            
            sumviszy = sumviszy + mspace(8,pairs(j,i)) * ( eirhoizy + &
                 epsilon(5,pairs(j,i))/rhoj2 ) * dwdx(2,j,i)
            
            sumviszz = sumviszz + mspace(8,pairs(j,i)) * ( eirhoizz + &
                 epsilon(6,pairs(j,i))/rhoj2 ) * dwdx(3,j,i)            
         enddo

         viscforce(1,i) = muc * mu * sumvisxx
         viscforce(2,i) = muc * mu * sumvisxy
         viscforce(3,i) = muc * mu * sumvisxz
         viscforce(4,i) = muc * mu * sumvisyx
         viscforce(5,i) = muc * mu * sumvisyy
         viscforce(6,i) = muc * mu * sumvisyz
         viscforce(7,i) = muc * mu * sumviszx
         viscforce(8,i) = muc * mu * sumviszy
         viscforce(9,i) = muc * mu * sumviszz
         
      enddo
      
    end subroutine viscous_force
