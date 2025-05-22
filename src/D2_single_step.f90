
      subroutine D2_single_step(mspace,ntype,itimestep,dt,av)
!c----------------------------------------------------------------------


!c Subroutine to determine the right hand side of a differential
!c equation in a single step for performing time integration

!c In this routine and its subroutines the SPH algorithms are performed.

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

!c     npairs    : maximum number of pairs interaction             [out]
!c     nfilas    : number of interaction for each fluid particle   [out]
!c     mrij      : matrix of ri-rj for all fluid particles for each interaction [out]
!c     mxiij      : matrix of xi-xj for all fluid particles for each interaction [out]         
!c     mvij      : matrix of vi-vj for all fluid particles for each interaction [in]
!c     mvxij     : matrix of vxi-vxj for all fluid particles for each interaction [in]      

!c     t         : Temperature                                      [in/out]
!c tdsdt     : Production of viscous entropy t*ds/dt               [out]
!c av        : Monaghan average velocity                           [out]

      
      
!c rdomain   : smoothing length                                     [in]
      
      implicit none
      include 'param.inc'

      integer i,j,itimestep, ntotal,npairs,ntype(2),nfluid,nvirt
!c     parameter ( npairs = (kappa0+2)**dim)
      parameter ( npairs = (kappa0+5)**dim)
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      double precision t(fintime+1),dt,rdomain,av
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      integer pairsv(npairs,ntype(1)),nfilasv(ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1)), &
           mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1)),          &
           mrijv(npairs,ntype(1)),mxijv(3,npairs,ntype(1))
      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))
      double precision rrr,xi,yi,dx(dim),r,dwdx0(dim),ww

!      double precision start_x, start_z, vel, length, angle
      integer num_points,nstepsf

      rdomain = kappa0 * h0 * 1.00000007 !equivalent to hsml

      do i=1,ntype(1)
         mspace(13,i) = h0      !rdomain
      enddo

      if (mod(itimestep,save_step).eq.0) then
         write(*,*)'-------------------------------------'
         write(*,*)'Max number of n pairs = ',npairs
         write(*,*)'kappa0 = ',kappa0,'  h0 = ',h0
         write(*,*)'rdomain = ',rdomain
         write(*,*)'-------------------------------------'
      endif
      

!      call  input(mspace,ntotal,nfluid,nvirt,1)


!      vel = 0.5
!      nstepsf = int(sqrt(2*h0/g)/dt0)
!      write(*,*)'nstepsf = ',nstepsf

!      write(*,*)'ntype',ntype
!      do i=1,ntype(1)-ntype(2)+1
!         write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!              ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!              ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!              ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
!      enddo      
      
!      if (mod(itimestep, nstepsf) == 0) then
!         start_x = h0*3
!         start_z = h0*3
!         length = 0.1
!         angle = 0.
!         num_points = INT(length / h0) + 1
!         write(*,*)'Number of new particles = ',num_points
!         call fuente(start_x, start_z, vel, angle, num_points, ntype, mspace)
!         write(*,*)'Fuente modifique ntype  ',ntype
!      endif

!      write(*,*)'ntype 1 ',ntype
!      do i=1,ntype(1)-ntype(2)+1
!         write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!              ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!              ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!              ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
!      enddo
         
!         stop

         
      if (mod(itimestep, 10) == 0) then
         write(*, '(A)', advance='no') '.'
      endif


!      mrij = 0.
!      mxij = 0.
!      mvij = 0.
!      mvxij = 0.

!      write(*,*)'11'
      
      if (cns == 1)then
      !   call neighboring_search(rdomain,mspace,ntype,npairs,pairs,nfilas, &
      !        mrij,mxij,mvij,mvxij)
      endif

      if (cns == 2)then
      !   call neighboring_search_P(rdomain,mspace,ntype,npairs,pairs,nfilas, &
      !        mrij,mxij,mvij,mvxij)
      endif

      if (cns == 3)then
         call D2_neighboring_Grid_Hashing(rdomain,mspace,ntype,npairs,pairs, &
              nfilas,mrij,mxij,mvij,mvxij,itimestep)
      endif

      if (cns == 4)then
      !   call neighboring_Grid_Hashing_P(rdomain,mspace,ntype,npairs,pairs, &
      !        nfilas,mrij,mxij,mvij,mvxij)
      endif

      if (cns == 5)then
      !   call neighboring_Grid_Hashing_FVC(rdomain,mspace,ntype,npairs,pairs, &
      !        nfilas,mrij,mxij,mvij,mvxij,itimestep)
      endif

!      write(*,*)'ntype neighbor ',ntype
!      do i=1,ntype(1)-ntype(2)+1
!         write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!              ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!              ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!              ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
!      enddo
      
!      write(*,*)'22'
      call wijdwij(mspace,pairs,mrij,mxij,npairs,nfilas,ntype,w,dwdx) !2D y 3D

!      write(*,*)'ntype wijdwij ',ntype
!      do i=1,ntype(1)-ntype(2)+1
!         write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!              ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!              ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!              ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
!      enddo
      
      
!      write(*,*)'33'
      call density(mspace,ntype,npairs,pairs,nfilas,w) !2D y 3D

!      write(*,*)'ntype density ',ntype
!      do i=1,ntype(1)-ntype(2)+1
!         write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!              ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!              ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!              ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
!      enddo
!c      write(*,*)'44'
!c      do i=1,nmax
!c         do j=1,nfilas(i)
!c            write(*,*)i,j,pairs(j,i),w(j,i),
!c     +           mspace(2,i),mspace(4,i),
!c     +           mspace(2,pairs(j,i)),mspace(4,pairs(j,i)),mrij(j,i)
!c         enddo
!c      enddo

!      write(*,*)'44'
      call presioni(mspace,ntype, itimestep) !2D y 3D

!      write(*,*)'ntype presioni ',ntype
!      do i=1,ntype(1)-ntype(2)+1
!         write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!              ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!              ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!              ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
!      enddo

!      write(*,*)'55'
      call momento(mspace,ntype,npairs,pairs,nfilas,w,dwdx,mrij,mxij, &
          mvij,mvxij,pairsv,nfilasv,mrijv,mxijv)

!      write(*,*)'ntype momento ',ntype
!      do i=1,ntype(1)-ntype(2)+1
!         write(*,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(4,i)) &    
!              ,real(mspace(5,i)),real(mspace(7,i)),real(mspace(8,i))  &
!              ,real(mspace(9,i)),real(mspace(10,i)),real(mspace(11,i)) &
!              ,int(mspace(12,i)),real(mspace(13,i)),real(mspace(26,i))
!      enddo
      
!c      do i=1,ntype(1)
!c         do j=1,nfilas(i)
!c            dx(1) = mspace(2,i)-mspace(2,pairs(j,i))
!c            dx(2) = mspace(4,i)-mspace(4,pairs(j,i))
!c            rrr = sqrt(dx(1)**2 + dx(2)**2)
!c            call kernel(rrr,dx,mspace(13,i),ww,dwdx0)
!c            if(i.eq.450)write(*,*)i,j,pairs(j,i),mspace(2,pairs(j,i)),
!c     +           mspace(4,pairs(j,i)),w(j,i),ww,dwdx(1,j,i),dwdx0(1),
!c     +           dwdx(3,j,i),dwdx0(2),rrr
!c     enddo
!c        write(*,*),mspace(20,i),mspace(21,i),mspace(22,i)
!c      enddo
      

      
    end subroutine D2_single_step


