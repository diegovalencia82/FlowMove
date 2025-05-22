
ulimit -s unlimited

export OMP_DISPLAY_ENV=TRUE

export OMP_DYNAMIC=TRUE
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

export OMP_STACKSIZE=8G

export OMP_WAIT_POLICY=ACTIVE

export OMP_NUM_THREADS=2



../sph
