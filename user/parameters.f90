module parameters

    use constants
    implicit none

#ifdef MPIP
    include "mpif.h"
#endif
                                              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter  :: LB = BC_REFLECTIVE    !!                           Boundaries                           !!
    integer, parameter  :: RB = BC_REFLECTIVE    !!  LB = Left Boundary   DB = Down Boundary  FB = Front Boundary  !!
    integer, parameter  :: TB = BC_REFLECTIVE    !!  RB = Right Boundary  TB = Top Boundary   BB = Back Boundary   !!
    integer, parameter  :: DB = BC_REFLECTIVE    !!           (x)                 (y)                  (z)         !!
    integer, parameter  :: FB = BC_REFLECTIVE    !!                        Type of boundary                        !!
    integer, parameter  :: BB = BC_REFLECTIVE    !!  1 = Outflow    2 = Reflective    3 = Periodic    4 = Inflow   !!
                                              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer, parameter  :: limtr = LIMITER_VANLEER  !! Limiter

    integer, parameter  :: ndim = 3           !! Number of dimensions of the problem (Currently only supporting 3).
    integer, parameter  :: neq = 5            !! Number of equations to solve (Currently only supporting 5).
    integer, parameter  :: size = 6           !! 
    integer, parameter  :: nxtot = 50!4*2**size  !! Number of cells in each axis. If using self-gravity, nx, ny and nz
    integer, parameter  :: nytot = 50!4*2**size  !! MUST be a power of 2 with `size` it's power: nx, ny, nz = 2**size.
    integer, parameter  :: nztot = 50!4*2**size  !! If using MPI and self-gravity, nx, ny and nz MUST be a power of 2
    integer, parameter  :: choice = 2         !! like so: nx, ny, nz = (mpix, mpiy, mpiz)*2**size

    integer, parameter      :: num_out = 10  !! Number of output files 
    integer, parameter      :: outfile = VTK  !! Type of output. Options are VTK, DAT, BIN
    character(*), parameter :: outputpath = './DATA/'
#ifdef COOL
    character(*), parameter :: cooling_file = '../Z1.0.dat' !! Cooling table file location for atomic cooling
#endif
    character(*), parameter :: posfile = '/storage5/luis.arcos/25SG/posest75.dat' !! Star positions file for winds and point_gravity modules
    integer, parameter      :: nghost = 2   !! Order

    real, parameter  :: xl = -1.0!*PC           !! Position of first physical cell in the x axis
    real, parameter  :: xr = 1.0!*PC            !! Position of last physical cell in the x axis
    real, parameter  :: yl = -1.0!*PC           !! Position of first physical cell in the y axis
    real, parameter  :: yr = 1.0!*PC            !! Position of last physical cell in the y axis
    real, parameter  :: zl = -1.0!*PC           !! Position of first physical cell in the z axis
    real, parameter  :: zr = 1.0!*PC            !! Position of last physical cell in the z axis
    real, parameter  :: gamma = 5.0/3.0        !! Specific heat ratio. C_p/C_v
    real, parameter  :: mu0 = 1.0              !! Mean particule mass
    real, parameter  :: Gconst = GR            !! Gravitational constant. If using physical units, GR = 6.67430e-8 dyn cm^2/g^2
    real, parameter  :: cfl = 0.4              !! Courant-Friedrichs-Lewy condition. 0 < cfl < 1
    real, parameter  :: eta = 0.5e-2           !! Artificial Viscosity. 0 < eta < 0.5
    real, parameter  :: tfin = 1!1.9e5*YEAR      !! Total simulated time.
    real, parameter  :: cooling_limit = 0.5
                !!! MPI Constants definitions !!!
#ifdef MPIP          
    integer, parameter  :: mpix = 4            !! Number of subdivisions of the computational domain in the x axis
    integer, parameter  :: mpiy = 4            !! Number of subdivisions of the computational domain in the y axis
    integer, parameter  :: mpiz = 4            !! Number of subdivisions of the computational domain in the z axis 
    integer, parameter  :: np = mpix*mpiy*mpiz !! Total number of processor cores needed for the given subdivisions
    integer, parameter  :: nx = nxtot/mpix     !! Number of cells in each core in the x axis
    integer, parameter  :: ny = nytot/mpiy     !! Number of cells in each core in the y axis
    integer, parameter  :: nz = nztot/mpiz     !! Number of cells in each core in the z axis
#else
    integer, parameter  :: mpix = 1
    integer, parameter  :: mpiy = 1
    integer, parameter  :: mpiz = 1
    integer, parameter  :: nx = nxtot
    integer, parameter  :: ny = nytot
    integer, parameter  :: nz = nztot
    integer, parameter  :: np = 1
#endif

#ifdef MPIP
#ifdef DOUBLEP
    integer, parameter :: mpi_real_kind=mpi_real8
#else
    integer, parameter :: mpi_real_kind=mpi_real4
#endif
#else
    integer, parameter :: mpi_real_kind=0
#endif

    real, parameter  :: dtout = tfin/num_out   !! Time step between outputs

    integer, parameter  :: nxmax = nx+nghost    !!
    integer, parameter  :: nymax = ny+nghost    !!
    integer, parameter  :: nzmax = nz+nghost    !! Total number of cells including
                                                !! ghost cells in each core given
    integer, parameter  :: nxmin = 1-nghost     !! the order and reconstruction scheme
    integer, parameter  :: nymin = 1-nghost     !!
    integer, parameter  :: nzmin = 1-nghost     !!

    real, parameter     :: dx = (xr-xl)/nxtot   !! cell size in the x axis
    real, parameter     :: dy = (yr-yl)/nytot   !! cell size in the y axis
    real, parameter     :: dz = (zr-zl)/nztot   !! cell size in the z axis
#ifdef GRAV
#ifdef MPIP
    integer, parameter  :: lvlm = size
#else
    integer, parameter  :: lvlm = size
#endif
#endif
    logical, parameter      :: logged = .false.            !! Sets whether or not create a logging file
    integer, parameter      :: logu = 12
    character(*), parameter :: logfile='./log125SGG.dat'   !! Name of the logging file

    character(*), parameter :: paramFile = "tlaloque.dat"  !! Name of the file to print all runtime parameters
end module parameters
