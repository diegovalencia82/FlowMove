
      subroutine tasa_deformacion_epsilon(mspace,ntype,npairs,pairs, &
          nfilas,w,dwdx,mxij,mvij,mvxij,epsilon)
      
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
      
      implicit none
      include 'param.inc'
        
      integer i,j,npairs,ntype(2)
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1)), &
           mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))
      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))
      double precision epsilon(6,ntype(1)),mjrhoj,sumexx,sumeyy,sumezz
      double precision sumexy1,sumexy2,sumexz1,sumexz2,sumeyz1,sumeyz2, &
           divvdwdx,sumdivg,c23
      
      c23 = 2.0d0 / 3.0d0
      
      do i=1,ntype(1)-ntype(2)
         sumexx = 0.0
         sumeyy = 0.0
         sumezz = 0.0
         sumexy1 = 0.0
         sumexy2 = 0.0
         sumexz1 = 0.0
         sumexz2 = 0.0
         sumeyz1 = 0.0
         sumeyz2 = 0.0
         sumdivg = 0.0
         do j=1,nfilas(i)
            mjrhoj = mspace(8,pairs(j,i))/mspace(9,pairs(j,i))
            divvdwdx = -( mvxij(1,j,i)*dwdx(1,j,i) + &
                mvxij(2,j,i)*dwdx(2,j,i) + mvxij(3,j,i)*dwdx(3,j,i) )
            sumdivg = sumdivg + mjrhoj*divvdwdx
            
            sumexx = sumexx + mjrhoj*(-mvxij(1,j,i))*dwdx(1,j,i)
            sumeyy = sumeyy + mjrhoj*(-mvxij(2,j,i))*dwdx(2,j,i)
            sumezz = sumezz + mjrhoj*(-mvxij(3,j,i))*dwdx(3,j,i)
            
            sumexy1 = sumexy1 + mjrhoj*(-mvxij(2,j,i))*dwdx(1,j,i)
            sumexy2 = sumexy2 + mjrhoj*(-mvxij(1,j,i))*dwdx(2,j,i)

            sumexz1 = sumexz1 + mjrhoj*(-mvxij(3,j,i))*dwdx(1,j,i)
            sumexz2 = sumexz2 + mjrhoj*(-mvxij(1,j,i))*dwdx(3,j,i)

            sumeyz1 = sumeyz1 + mjrhoj*(-mvxij(3,j,i))*dwdx(2,j,i)
            sumeyz2 = sumeyz2 + mjrhoj*(-mvxij(2,j,i))*dwdx(3,j,i)
         enddo
         epsilon(1,i) = 2.*sumexx - c23*sumdivg     !e(xx)
         epsilon(2,i) = sumexy1 + sumexy2           !e(xy) = e(yx)
         epsilon(3,i) = sumexz1 + sumexz2           !e(xz) = e(zx)
         epsilon(4,i) = 2.*sumeyy - c23*sumdivg     !e(yy)
         epsilon(5,i) = sumeyz1 + sumeyz2           !e(yz) = e(zy)
         epsilon(6,i) = 2.*sumezz - c23*sumdivg     !e(zz)

      enddo
      


    end subroutine tasa_deformacion_epsilon
