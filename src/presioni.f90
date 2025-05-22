
      subroutine presioni(mspace,ntype,itimestep)

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

      implicit none
      include 'param.inc'

      integer i, npairs, ntype(2), itimestep
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      double precision b, beta0

      
!c     Artificial equation of state for the artificial compressibility
!c     rho : Density                                                 [in]
!c     u   : Internal energy                                         [in]
!c     p   : Pressure                                               [out]
!c     c   : sound velocity                                         [out]
!c     b   : is a problem dependent parameter, which sets a limit
!c           for the maximum change of the density.
!c           In most circumstances, b can be taken as the initial
!c           pressure (Morris et al., 1997; Schlatter, 1999). 

!c     Equation of state for artificial compressibility

!c      implicit none
!c      double precision rho, u, p, c, b
!c      double precision gamma, rho0
!c      double precision ht,y,rho1,beta,g

!c      ht = 1.e-3*400
!c      g=9.8
!c      gamma=7.
!c      beta=5.0!0.003
!c      rho0=1000.
!c      c = 1480
      
      
!c     Artificial EOS, Form 1 (Monaghan, 1994)
!c     write(*,*)'Artificial EOS, Form 1 (Monaghan, 1994)'
!c      b = beta*g*ht*rho0/gamma !1.013e5
!c      p = b*((rho/rho0)**gamma-1)
!c      c = sqrt(beta*g*ht)!1480.

!c     Artificial EOS, Form 2 (Morris, 1997)
!c      write(*,*)'Artificial EOS, Form 2 (Morris, 1997)'
!c      c = 0.01
!c      p = c**2 * rho

!c     Artifical EOS, form3 (DualPhysics)
!c      beta = 10.
!c      c = beta*sqrt(9.8*ht)
!c      b = c*c*rho0/gamma
      !c      p = b*((rho/rho0)**gamma-1)

      if(rbeta.eq.1)then
         beta0=beta
      else
         beta0 = 1+itimestep*rbeta*5 
         if(beta0.gt.beta)beta0=beta         
      endif

      do i=1,ntype(1)-ntype(2)
         mspace(14,i) = beta0*sqrt(2*g*ht)
         b = rho0*mspace(14,i)*mspace(14,i)/gamma
         mspace(10,i) = b*((mspace(9,i)/rho0)**gamma-1.)
      enddo
      
    end subroutine presioni
