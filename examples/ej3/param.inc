!-------------------------------------------------------
!     Including file for parameters and constants used
!     in the entire SPH software packages.
!-------------------------------------------------------

      double precision pi
      parameter ( pi = 3.14159265358979323846 )

!     dim : Dimension of the problem (1, 2 or 3)
      integer dim
      parameter ( dim = 2)

!     nmax	: Maximum number of particles
      integer nmax
      parameter ( nmax = 100000 ) 
!

!-------------------------------------------------------
!     dominio de simulación
      double precision dxmin, dxmax, dymin, dymax, dzmin, dzmax
      parameter ( dxmin = -0.100000, dxmax = 1.710000 )
!      parameter ( dymin = -0.100000, dymax = 0.250000 )
      parameter ( dzmin = -0.100000, dzmax = 0.80000 )

!-------------------------------------------------------
!      call neighboring_search
!             cns = 1, direct Search for neighbors
!             cns = 2, direct parallelized neighbor search
!             cns = 3, Grid Hashing Search for neighbors
!             cns = 4, Grid Hashing parallelized Search for neighbors
!             cns = 5, Grid Hashing Search for neighbors, meausre raidus for aparent viscosity 
      integer cns
      parameter (cns = 3)

!    if cns = 4 you have to set the number of processors
     integer np
     parameter (np = 8)

!-------------------------------------------------------
      integer maxtimestep, initime, fintime, inioutfile
      parameter ( initime = 0, fintime = 200000 )  !This refers to the number of steps itself, not the time in seconds.
      parameter ( inioutfile = 000)	
      character infile*80, infilev*80
      parameter ( infile = 'snapshot_000' )
      parameter ( infilev = 'snapshot_000' ) ! carga el archivo cuando inicia con archivos sin particulas de frontera (Siempre carga el archivo inicial el que tiene todas las particulas

!-------------------------------------------------------
!     Time step integration
      double precision dt0, pdt0
      parameter ( dt0 = 0.0001, pdt0 = 0.73 )

!-------------------------------------------------------
!      To save steps
       integer save_step
       parameter ( save_step = 100)

!-------------------------------------------------------
!     Definition of the support domain
!     h0 is taken from the space of particles in the initial conditions
      double precision h0, kappa0
      parameter ( kappa0 = 2.0d0, h0 = 0.006d0 ) ! h0 = 0.01

!-------------------------------------------------------
!     Smoothing kernel function
!     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!         = 2, Gauss kernel (Gingold and Monaghan 1981)
!         = 3, Quintic kernel (Morris 1997)
      integer skf
      parameter(skf = 1)

!-------------------------------------------------------
!     ht 	    Maximum high of the block, it is equal to the initial conditions
!     The measurement must be in meters
      double precision ht
      parameter ( ht = 0.3d0)

!-------------------------------------------------------
!     g = 9.82. Gravity
      double precision g
      parameter ( g = 9.82d0 )

!-------------------------------------------------------
!     Artificial EOS.
!     Initial density of particles, rho0
!     rbeta es un parametro de regularizacion
      double precision gamma, beta, rho0, rbeta
      parameter ( gamma=7.0, beta = 1.0, rbeta = 1., rho0 = 1000.0)
!      parameter ( gamma=7.0, beta = 2.23607, rbeta = 0.05, rho0 = 1000.0)

!-------------------------------------------------------
!     viscous force
      double precision mu, muc
      parameter ( mu=1.002e-3, muc=100)
!      parameter ( mu=1.002e-3, muc=300)		
!      parameter ( mu=8.41, muc=1.)
!      parameter ( mu=800.41, muc=1.)	


!-------------------------------------------------------
!      mvisc = 1 constant viscosity
!      mvisc = 2 non-constant viscosity
!      nn = índice de comportamiento reológico o índice de fluidez.
!      mvisK = 1 indica parametro de viscosidad K constante en el modelo de Herschel-Bulkley
!      mvisK = 2 indica parametro de viscosidad K depende de la altura en el modelo de Herschel-Bulkley
      integer mvisc, mvisk
      double precision nn
      parameter (mvisc = 1, nn = 1.0, mvisk = 1)

!     if mvisk = 2 you have to choose
!     phi_base, phi_surface
      double precision phi_base, phi_surface, intrinsic_viscosity
      parameter( phi_base = 0.8d0, phi_surface = 0.5d0, intrinsic_viscosity = 2.5d0)

!-------------------------------------------------------

!      source_mode = y:  The code use the source mode
!      mu_new = viscosity of soruce mode fluid
!      QQ = flow rate value [m^2 / s]
!      length = length of the line source of fluid
!      start_x = initial position of the source in x axe
!      start_z = initial position of the source in y axe
!      angle = source angle
      character source_mode
      parameter (source_mode = 'y')
      double precision QQ, mu_new, start_x, start_y, start_z, angle, v_fuente0, d_horizontal
      parameter (QQ = 0.03d0, mu_new = 1.002e-3)	
      parameter (start_x = 0.02, start_y = 0.2, start_z = 0.1, angle = 90.0d0)	
      parameter (v_fuente0 = 1.0d0, d_horizontal = 0.2d0)     	   

