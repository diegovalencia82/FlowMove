
      subroutine viscous_model(mspace,ntype,epsilon2,mud)

!c     mspace(10,i) = p-- pressure of particles                     [in]
!c     mspace(26,i) = eta_c     : Coeficiente de Viscosidad No Lineal [in]      
      
      implicit none
      include 'param.inc'

      integer i,ntype(2)
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      double precision epsilon2(ntype(1)),gammap
      double precision c,phi,nc,mud(ntype(1)),nume,dem,K,tauB,mu_inf, &
           mu00,nn1, eta_m,aux!, concentracion
      double precision concentration_profile, krieger_dougherty, maron_pierce, htt

      
!      c = 10.
!      phi = 30 * pi / 180.
!      mu_inf = mu
!      mu00 = 50.*mu_inf
      
!      do i=1,ntype(1)
!         mud(i) = mu * 50
!      enddo

      ! Modelo de Binghman
!c      do i = 1,ntype(1)!-ntype(2)
!c         gammap = sqrt(0.5*epsilon2(i))
!c         tauB = ( c + mspace(10,i)*tan(phi) ) * 0.5
!c         K = mu00 / tauB
!c         nume = K*mu_inf*gammap + mu00
!c         dem = K*gammap + 1
!c         mud(i) = nume / dem
!c         write(*,*)i,mspace(12,i),mspace(10,i),K,mu_inf,gammap,mu00,
!c     +        mud(i),tauB
!c      enddo

!c      nc = mu / 1000.
!c      concentracion = 90 !%
!c     p = (-1./100.)*concentracion + 1
!c      p = 1.!1./100000.
!c      nc = mu * p


      
      if (mvisk.eq.1)then
         do i = 1,ntype(1)-ntype(2)
            !nn1 = nn
            nn1 = mspace(19,i)
            mud(i) = mspace(26, i) * (0.5*epsilon2(i))**((nn1-1)/2)
         enddo
      endif

      if (mvisk.eq.2)then
         htt = h0 * 2
         !aux = 0
         do i = 1,ntype(1)-ntype(2)
            !concentracion = concentration_profile(mspace(17,i), htt, phi_base, phi_surface)
            mspace(18,i) = concentration_profile(mspace(17,i), htt, phi_base, phi_surface)
            !eta_m = krieger_dougherty(concentracion, phi_base, mspace(26,i), intrinsic_viscosity)
            eta_m = krieger_dougherty(mspace(18,i), phi_base, mspace(26,i), intrinsic_viscosity)
            !eta_m = maron_pierce(phi(i), phi_max, mspace(26, i))
            if(eta_m.lt.mspace(26,i))eta_m=mspace(26,i)
            !nn1 = nn
            nn1 = mspace(19,i)
            mud(i) = eta_m * (0.5*epsilon2(i))**((nn1-1)/2)
            !write(*,*)i,eta_m,mspace(18,i),mspace(17,i)
            !if(aux.lt.eta_m)aux = eta_m
         enddo
         do i = ntype(1)-ntype(2)+1,ntype(1)
            mud(i) = mspace(26,i)
         enddo
      endif

      !write(*,*)'eeeeeeeeeeeee',aux

      
    end subroutine viscous_model


    ! =================================================================================


        ! 1. Perfil lineal de concentración de sólidos
    function concentration_profile(height, max_height, phi_base, phi_surface) result(phi)
        double precision, intent(in) :: height, max_height, phi_base, phi_surface
        double precision :: phi
        
        phi = phi_base - (phi_base - phi_surface) * (height / max_height)
    end function concentration_profile

    ! 2. Modelo de Krieger-Dougherty
    function krieger_dougherty(phi, phi_max, eta_0, intrinsic_viscosity) result(eta)
        double precision, intent(in) :: phi, phi_max, eta_0, intrinsic_viscosity
        double precision :: eta, ratio
        
        ! Evitar división por cero o valores negativos
        ratio = max(1.0d-12, 1.0d0 - phi / phi_max)
        eta = eta_0 * ratio**(-intrinsic_viscosity * phi_max)
    end function krieger_dougherty

    ! 3. Modelo de Maron-Pierce
    function maron_pierce(phi, phi_max, eta_0) result(eta)
        double precision, intent(in) :: phi, phi_max, eta_0
        double precision :: eta, ratio
        
        ! Evitar división por cero o valores negativos
        ratio = max(1.0d-12, 1.0d0 - phi / phi_max)
        eta = eta_0 * ratio**(-2.0d0)
    end function maron_pierce
