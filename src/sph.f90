!module shared_data
!    implicit none
!    double precision, allocatable :: mspace(:,:)
!end module shared_data
  
      program SPH

!c------------------------------------------------------------------
!c     This is a three dimensional SPH code. the followings are the 
!c     basic parameters needed in this codeor calculated by this code

!c     ntotal-- total particle number ues                             [in]
      
!c     dt--- Time step used in the time integration                   [in]
!c     time -- time of the snapshot
!c     itimestep -- number of the time step.
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

      
!      use shared_data
      implicit none
      include 'param.inc'
      
      integer ntotal,nfluid,nvirt,ntype(2) !, itype(nmax),id(nmax)
      double precision dt
!      double precision mspace(26,nmax)
      !      double precision, dimension(26, nmax) :: mspace
      double precision, allocatable :: mspace(:,:)
      
      integer :: start_clock, end_clock, rate
      real :: total_time

      ! Asignar memoria
      allocate(mspace(26, nmax))
      
! Obtener el tiempo inicial
   
      call SYSTEM_CLOCK(start_clock, rate)

      dt = dt0

      write(*,*)'********** Start ***********'
      call input(mspace,ntotal,nfluid,nvirt,1)
      call input(mspace,ntotal,nfluid,nvirt,0)
      write(*,*)'********** Start ***********'
      
      ntype(1) = nfluid + nvirt
      ntype(2) = nvirt
      write(*,*)' **************************************************'
      write(*,*)'        The  maximal time steps = ', fintime+1
      write(*,*)'        Time steps integration  = ', dt
      write(*,*)' **************************************************' 

      if(dim == 2)call D2_time_integration(mspace,ntotal,ntype,dt)
      if(dim == 3)call time_integration(mspace,ntotal,ntype,dt)


      deallocate(mspace)  ! Libera memoria al final
      
! Obtener el tiempo final
      call SYSTEM_CLOCK(end_clock, rate)
! Calcular el tiempo de CPU consumido
      total_time = real(end_clock - start_clock) / real(rate)
      open(83,file='time_CPU.dat')
      write(83,*) 'START TIME (clock ticks) = ', start_clock
      write(83,*) 'END TIME (clock ticks) = ', end_clock
      write(83,*) 'Total time (minutes) = ', total_time / 60
      close(83)
      
      end program SPH
