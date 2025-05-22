
      subroutine input(mspace,ntotal,nfluid,nvirt,sel)

!c-----------------------------------------------------------------

!c     time -- time of the snapshot
!c     itimestep -- number of the time step.
!c     mspace(1,i) = id -- label of each particle                   [out]      
!c     mspace(2 to 4,i) = x-- coordinates of particles              [out]
!c     mspace(5 to 7,i) = vx-- velocities of particles              [out]
!c     mspace(8,i) = mass-- mass of particles                       [out]
!c     mspace(9,i) = rho-- dnesities of particles                   [out]
!c     mspace(10,i) = p-- pressure of particles                     [out]
!c     mspace(11,i) = u-- internal energy of particles              [out]
!c     mspace(12,i) = itype-- types of particles                    [out]
!c     mspace(13,i) = hsml-- smoothing lengths of particles         [out]
!c     mspace(14,i) = c-- sound velocity of particle                [out]
!c     mspace(15,i) = s-- entropy of particles                      [out]
!c     mspace(16,i) = e-- total energy of particles                 [out]      
!c     mspace(17 ,i) = Distancia i fluido-frontera                  [out]
!c     mspace(18 ,i) = Viscosidad K dependiente de la concentra     [out]
!c     mspace(19 ,i) = Exponente n para Hercshel-Bulkley            [out]
!c     mspace(20 to 22,i) = dvx = dvx/dt, force per unit mass       [out]
!c     mspace(23,i) = du        : du = du/dt                        [out]
!c     mspace(24,i) = ds        : ds = ds/dt                        [out]
!c     mspace(25,i) = drho      : drho = drh,o/dt                   [out]
!c     mspace(26,i) = eta_c     : drho = drh,o/dt                   [out]      Coeficiente de Viscosidad No Lineal
!c     ntotal-- total particle number                               [out]

      implicit none
      include 'param.inc'
      
      integer itype(nmax),id(nmax)
      integer ntotal,nfluid,nvirt,itimestep
      double precision time,x(dim, nmax), vx(dim, nmax), mass(nmax),&
           p(nmax), u(nmax), hsml(nmax), rho(nmax),etac(nmax),n_exp(nmax),&
           visck(nmax)
      !double precision mspace(26,nmax)
      double precision, dimension(26, nmax) :: mspace
      integer i, i1, i2, d, im, sel
      
!c     load initial particle information from external disk file
      !c     open(10,file="snapshot_000",status='old')

      if(sel.eq.0)open(10,file=infile,status='old')
      if(sel.eq.1)open(10,file=infilev,status='old')
      

      if(sel.eq.0)write(*,*)' ****************************************'
      if(sel.eq.0)write(*,*)'Loading initial particle configuration...'
      read (10,*) itimestep,time,ntotal,nfluid,nvirt

      i=1
      if(sel.eq.0)then
         write(*,*)'       Total number of particles       ', ntotal
         write(*,*)'       Total number of fluid particles       ', nfluid
         write(*,*)'       Total number of virtual particles       ', nvirt
         write(*,*)'       itimestep = ',itimestep,'   Time = ',time
         write(*,*)'   ************************************************'
      endif

 11   read(10,*,end=12)id(i), (x(d, i),d = 1, dim), &
           (vx(d, i),d = 1, dim),mass(i), rho(i), p(i), u(i),itype(i), &
           hsml(i), etac(i), viscK(i) ,n_exp(i)
      i=i+1
      goto 11
 12   close(10)
      i=i-1


      
      if(sel.eq.0)then
         i1 = 1
         i2 = nfluid !ntotal
      endif
      if(sel.eq.1)then
         i1 = nfluid + 1
         i2 = ntotal
      endif

      do i=i1,i2
         if(dim.eq.2)then 
            mspace(1,i)  = id(i)
            mspace(2,i)  = x(1,i)
!c     mspace(3,i)  = x(2,i)
            mspace(4,i)  = x(2,i)
            mspace(5,i)  = vx(1,i)
!c     mspace(6,i)  = vx(2,i)
!            if(sel.eq.0)then
               mspace(7,i)  = vx(2,i)
               mspace(8,i)  = mass(i)
               mspace(9,i)  = rho(i)
               mspace(10,i)  = p(i)
               mspace(11,i)  = u(i)
               mspace(12,i) = itype(i)
               mspace(13,i) = hsml(i)
               mspace(26,i) = etac(i)
               mspace(19,i) = n_exp(i)
!            endif
         endif
         if(dim.eq.3)then
            mspace(1,i)  = id(i)
            mspace(2,i)  = x(1,i)
            mspace(3,i)  = x(2,i)
            mspace(4,i)  = x(3,i)
            mspace(5,i)  = vx(1,i)
            mspace(6,i)  = vx(2,i)
            mspace(7,i)  = vx(3,i)
!            if(sel.eq.0)then
               mspace(8,i)  = mass(i)
               mspace(9,i)  = rho(i)
               mspace(10,i)  = p(i)
               mspace(11,i)  = u(i)
               mspace(12,i) = itype(i)
               mspace(13,i) = hsml(i)
               mspace(26,i) = etac(i)
               mspace(19,i) = n_exp(i)
!            endif
         endif
        
      enddo
      
    end subroutine input
