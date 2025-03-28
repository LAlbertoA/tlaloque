  subroutine initmain(nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq)
    
    use globals
    use constants
    use parameters, only: np, nx, ny, nz, ndim, logged, logu, logfile
#ifdef GRAV
    use parameters, only: lvlm
#endif
#ifdef MPIP
    use parameters, only: mpix, mpiy, mpiz, LB, RB, DB, TB, FB, BB
#endif
    implicit none
    
#ifdef MPIP
    include "mpif.h"
#endif
    
    integer, intent(in)      :: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq
#ifdef GRAV
    integer                  :: nxl, nyl, nzl, level
#endif
    
    allocate(    dims(0:ndim-1))
    allocate(  coords(0:ndim-1))
    allocate( periods(0:ndim-1))

#ifdef MPIP
    call mpi_init(err)
    call mpi_comm_rank(mpi_comm_world,rank,err)
    call mpi_comm_size(mpi_comm_world,nprocs,err)
    
    if (nprocs.ne.np) then
       print*, 'processor number (',nprocs,') is not equal to pre-defined number (',np,')'
       call mpi_finalize(err)
       stop
    end if

    dims(0) = mpix
    dims(1) = mpiy
    dims(2) = mpiz

    periods(0) = .false.
    periods(1) = .false.
    periods(2) = .false.

    if (LB==BC_PERIODIC.or.RB==BC_PERIODIC) then
       periods(0) = .true.
    endif
    if (DB==BC_PERIODIC.or.TB==BC_PERIODIC) then
       periods(1) = .true.
    endif
    if (FB==BC_PERIODIC.or.BB==BC_PERIODIC) then
       periods(2) = .true.
    endif
    
    call mpi_cart_create(mpi_comm_world, ndim, dims, periods, 1, comm3d, err)
    call mpi_comm_rank(comm3d, rank, err)
    call mpi_cart_coords(comm3d, rank, ndim, coords, err)

    print '(A,i3,A,3(i3),A)', "Processor ", rank, " with coordinates ", coords(:), " active."

    call mpi_cart_shift(comm3d, 0, 1, left , right, err)
    call mpi_cart_shift(comm3d, 1, 1, down , top  , err)
    call mpi_cart_shift(comm3d, 2, 1, front, back , err)

    call mpi_barrier(comm3d, err)

    print '(A,i3,A,3i3,A,2i3,A,2i3,A,2i3,A)', "Hi, I'm rank ", rank, " with coordinates ", &
         coords(:), " and neighbors ", left, right, " in x, ", down, top, " in y, ", front, back, " in z"
#else
    dims(:) = 1

    periods(:) = .false.

    coords(:) = 0
    
    rank = 0
#endif
    allocate(    U(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
    allocate(   UP(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
    allocate(  UPP(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
    allocate(UPrim(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
    allocate(    F(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
    allocate(    G(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
    allocate(    H(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
#ifdef GRAV
    allocate(          S(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax))
    allocate( PHI(nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1))
    allocate(PHIP(nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1))
    allocate(PHIT(nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1))
    allocate(  Error(0:lvlm))
    allocate(Residue(0:lvlm))
    allocate(    ret(0:lvlm))
    do level=0,lvlm
       
       nxl = Int(nx/2**(level))
       nyl = Int(ny/2**(level))
       nzl = Int(nz/2**(level))
       allocate(Error(level)%data(0:nxl+1, 0:nyl+1, 0:nzl+1))
       allocate(Residue(level)%data(0:nxl+1, 0:nyl+1, 0:nzl+1))
       allocate(ret(level)%data(0:nxl+1, 0:nyl+1, 0:nzl+1))
       Error(level)%nxl = nxl
       Error(level)%nyl = nyl
       Error(level)%nzl = nzl
       Error(level)%data = 0.0
       ret(level)%nxl = nxl
       ret(level)%nyl = nyl
       ret(level)%nzl = nzl
       ret(level)%data = 0.0
       Residue(level)%nxl = nxl
       Residue(level)%nyl = nyl
       Residue(level)%nzl = nzl
       Residue(level)%data = 0.0
       
    end do
#endif
    if (logged .eqv. .true. .and. rank == 0) then
       open(unit=logu,file=logfile,status='unknown')
       write(logu,'(a)') 'LOGGING FILE'
    endif

    if (rank == 0) then
       call print_parameters()
    endif
    
  end subroutine initmain
  
  subroutine initflow()

    use globals, only: U
    use user
    
    implicit none

    call user_initconds(U)
    
  end subroutine initflow

  subroutine print_parameters()

    use parameters

    implicit none

    character(len=20)           :: LeftB, RightB, TopB, DownB, FrontB, BackB
    character(len=100)          :: buff
    
    open(74, file=paramFile, status="unknown")

    write(74,*) "************************************************************"
    write(74,*) "***************** PARAMETERS FROM TLALOQUE *****************"
    write(74,*) "************************************************************"
    write(74,*) ""
    write(74,*) " BOUNDARIES: "
    write(74,*) ""
    if (LB == 1) then
       LeftB = "Outflow"
    elseif (LB == 2) then
       LeftB = "Reflective"
    elseif (LB == 3) then
       LeftB = "Periodic"
    elseif (LB == 4) then
       LeftB = "Inflow"
    endif

    if (RB == 1) then
       RightB = "Outflow"
    elseif (RB == 2) then
       RightB = "Reflective"
    elseif (RB == 3) then
       RightB = "Periodic"
    elseif (RB == 4) then
       RightB = "Inflow"
    endif

    if (TB == 1) then
       TopB = "Outflow"
    elseif (TB == 2) then
       TopB = "Reflective"
    elseif (TB == 3) then
       TopB = "Periodic"
    elseif (TB == 4) then
       TopB = "Inflow"
    endif

    if (DB == 1) then
       DownB = "Outflow"
    elseif (DB == 2) then
       DownB = "Reflective"
    elseif (DB == 3) then
       DownB = "Periodic"
    elseif (DB == 4) then
       DownB = "Inflow"
    endif

    if (FB == 1) then
       FrontB = "Outflow"
    elseif (FB == 2) then
       FrontB = "Reflective"
    elseif (FB == 3) then
       FrontB = "Periodic"
    elseif (FB == 4) then
       FrontB = "Inflow"
    endif

    if (BB == 1) then
       BackB = "Outflow"
    elseif (BB == 2) then
       BackB = "Reflective"
    elseif (BB == 3) then
       BackB = "Periodic"
    elseif (BB == 4) then
       BackB = "Inflow"
    endif
    
    write(74,*) "Left Boundary: ", trim(LeftB), "Right Boundary: ", trim(RightB)
    write(74,*) "Down Boundary: ", trim(DownB), "Top Boundary: ", trim(TopB)
    write(74,*) "Front Boundary: ", trim(FrontB), "Back Boundary: ", trim(BackB)

    write(74,*) ""

    write(74,*) "************************************************************"
    write(74,*) "************************************************************"
    write(74,*) ""
    write(74,*) " GENERAL PARAMETERS: "
    write(74,*) ""
    if (limtr == 1) then
       buff = "None"
    elseif(limtr == 2) then
       buff = "Van Leer"
    elseif(limtr == 3) then
       buff = "Minmod"
    elseif(limtr == 4) then
       buff = "UMIST"
    elseif(limtr == 5) then
       buff = "Woodward"
    elseif(limtr == 6) then
       buff = "Superbee"
    end if

    write(74,*) "Limiter: ", trim(buff)
    write(74,*) "Courant-Friedrichs-Lewy condition (cfl): ", cfl
    write(74,*) "Artificial viscosity constant (eta): ", eta
    write(74,*) "Order of solver: ", choice
#ifdef MPIP
    write(74,*) "MPI usage: true. Number of proccessors used: ", np
#else
    write(74,*) "MPI usage: false"
#endif
#ifdef GRAV
    write(74,*) "Gravity: true"
#else
    write(74,*) "Gravity: false"
#endif
    write(74,*) "Number of output files: ", num_out
    if (outfile == 1) then
       buff = "DAT"
    elseif (outfile == 2) then
       buff = "VTK"
    elseif (outfile == 3) then
       buff = "BIN"
    end if
    write(74,*) "Type of output file: ", trim(buff)
    write(74,*) "Output files folder: ", outputpath
    write(74,*) ""
    write(74,*) "************************************************************"
    write(74,*) "************************************************************"
    write(74,*) ""
    write(74,*) " MESH PARAMETERS: "
    write(74,*) ""
    write(74,*) "Number of total physical cells in x dimension: ", nxtot
    write(74,*) "Number of total physical cells in y dimension: ", nytot
    write(74,*) "Number of total physical cells in z dimension: ", nztot

#ifdef MPIP
    write(74,*) "Number of physical cells in each processor in x dimension: ", nx
    write(74,*) "Number of physical cells in each processor in y dimension: ", ny
    write(74,*) "Number of physical cells in each processor in z dimension: ", nz
#endif
    write(74,*) ""

    write(74,*) "Position of first physical cell in the x axis: ", xl
    write(74,*) "Position of last physical cell in the x axis: ", xr
    write(74,*) "Position of first physical cell in the y axis: ", yl
    write(74,*) "Position of last physical cell in the y axis: ", yr
    write(74,*) "Position of first physical cell in the z axis: ", zl
    write(74,*) "Position of last physical cell in the z axis: ", zr

    write(74,*) ""

    write(74,*) "Cell size in the x axis: ", dx
    write(74,*) "Cell size in the y axis: ", dy
    write(74,*) "Cell size in the z axis: ", dz

    write(74,*) ""

    write(74,*) "************************************************************"
    write(74,*) "************************************************************"
    write(74,*) ""
    write(74,*) " PHYSICS PARAMETERS: "
    write(74,*) ""
    write(74,*) "Heat capacity ratio gamma: ", gamma
    write(74,*) "Mean particle mass mu: ", mu0
#ifdef COOL
    write(74,*) "Cooling: true"
    write(74,*) "Cooling file: ", cooling_file
#else
    write(74,*) "Cooling: false"
#endif
    write(74,*) "Total simulated time in seconds: ", tfin
#ifdef GRAV
    write(74,*) "Gravitational constant: ", Gconst
#endif
    close(74)
    
  end subroutine print_parameters
  
