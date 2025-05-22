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
  !subroutine fuente(start_x, start_z, vel, angle, num_points, ntype, mspace)
subroutine D2_fuente(vel, masa_particula, num_points, ntype, mspace, posy)
  implicit none
  include 'param.inc'

  ! Variables de entrada
!  double precision, INTENT(IN) :: start_x, start_z, angle
  INTEGER, INTENT(IN) :: num_points
  ! Variables de salida
  double precision, dimension(26, nmax) :: mspace
  !  double precision, dimension(26, nmax) :: mspaca
  double precision, allocatable :: mspaca(:,:)
  double precision, DIMENSION(num_points) :: x_coords, y_coords, z_coords
  ! Variables locales
  double precision angle_rad, vel, masa_particula, c, b, posy
  INTEGER i,j, ntype(2),nfluid,nfluidt

  ! Asignar memoria
  allocate(mspaca(26, nmax))  ! local matrix
  
  ! Convertir ángulo a radianes
  angle_rad = angle * 3.14159265358979323846d0 / 180.0d0

  !constantes fisicas
  c = beta * sqrt(2 * g * ht) !sqrt(beta * g * ht)!
  b = rho0 * c *c / gamma

  mspaca = mspace

  ! Generar puntos
  DO i = 1, num_points
     x_coords(i) = start_x + (i - 1) * 2.*h0 * COS(angle_rad)
     z_coords(i) = start_z + (i - 1) * 2.*h0 * SIN(angle_rad)
!     y_coords(i) = posy
  END DO

  nfluid = ntype(1)-ntype(2)
  nfluidt = nfluid + num_points
  j = 0
  do i = nfluid + 1, nfluidt
     j = j + 1
     mspace(1,i) = i
     mspace(2,i) = x_coords(j)
     mspace(3,i) = 0.!y_coords(j)
     mspace(4,i) = z_coords(j)
     mspace(5,i) = vel * sin(angle_rad)
     mspace(6,i) = 0.0d0
     mspace(7,i) = vel * cos(angle_rad)
     mspace(9,i) = rho0 !* (1 + (rho0 * g * (ht - abs(mspace(4,i) )) / b)) ** (1.0 / gamma)
     mspace(8,i) = masa_particula!mspace(9,i) * h0**3
     mspace(10,i) = b*((mspace(9,i)/rho0)**gamma - 1.)  
     mspace(11,i) = 357.1
     mspace(12,i) = 1
     mspace(13,i) = h0
     mspace(14,i) = c
!     mspace(15,i) =
!     mspace(16,i) =
!     mspace(17,i) =
!     mspace(18,i) =
!     mspace(19,i) =
     mspace(20,i) = mspaca(20,j)
     mspace(21,i) = mspaca(21,j)
     mspace(22,i) = mspaca(22,j)
!     mspace(23,i) =
!     mspace(24,i) =
!     mspace(25,i) =
     mspace(26,i) = mu_new
     write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
          ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
          ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
          ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
  enddo

  j = nfluid
  ntype(1) = ntype(1) + num_points

  do i = nfluidt + 1, ntype(1)
     j = j + 1
     mspace(1,i) = i
     mspace(2,i) = mspaca(2,j)
     mspace(3,i) = mspaca(3,j)
     mspace(4,i) = mspaca(4,j)
     mspace(5,i) = mspaca(5,j)
     mspace(6,i) = mspaca(6,j)
     mspace(7,i) = mspaca(7,j)
     mspace(8,i) = mspaca(8,j)
     mspace(9,i) = mspaca(9,j)
     mspace(10,i) = mspaca(10,j)
     mspace(11,i) = mspaca(11,j)
     mspace(12,i) = mspaca(12,j)
     mspace(13,i) = mspaca(13,j)
     mspace(14,i) = mspaca(14,j)
!     mspace(15,i) =
!     mspace(16,i) =
!     mspace(17,i) =
!     mspace(18,i) =
!     mspace(19,i) =
!     mspace(20,i) =
!     mspace(21,i) =
!     mspace(22,i) =
!     mspace(23,i) =
!     mspace(24,i) =
!     mspace(25,i) =
     mspace(26,i) = mspaca(26,j)
!     write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!          ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!          ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!          ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
  enddo


  deallocate(mspaca)  ! Libera memoria al final
END SUBROUTINE D2_FUENTE
  
