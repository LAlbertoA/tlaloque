module globals
  implicit none

  !!--------------BIG ARRAYS FOR HYDRO SOLVER--------------!!
  real, dimension(:,:,:,:), allocatable               :: U
  real, dimension(:,:,:,:), allocatable               :: UP
  real, dimension(:,:,:,:), allocatable               :: UPP
  real, dimension(:,:,:,:), allocatable               :: UPrim
  real, dimension(:,:,:,:), allocatable               :: F
  real, dimension(:,:,:,:), allocatable               :: G
  real, dimension(:,:,:,:), allocatable               :: H
  !!--------------PARAMETERS FOR HYDRO SOLVER--------------!!
  real                                                :: t, dt
  real                                                :: tout, dth
  integer                                             :: itprint, nout, it, png
#ifdef GRAV
  !!-------------BIG ARRAYS FOR POISSON SOLVER-------------!!
  real, dimension(:,:,:,:), allocatable               :: S
  real, dimension(:,:,:), allocatable                 :: PHI
  real, dimension(:,:,:), allocatable                 :: PHIP
  real, dimension(:,:,:), allocatable                 :: PHIT
  type :: MultiGridLevelData
     integer                                          :: nxl
     integer                                          :: nyl
     integer                                          :: nzl
     real*8, dimension(:,:,:), allocatable            :: data
  end type MultiGridLevelData
  type(MultiGridLevelData), dimension(:), allocatable :: Error, Residue, ret
#endif
  !!-----------------PARAMETERS FOR MPI-------------------!!
  integer, dimension(:), allocatable                  :: dims, coords
  logical, dimension(:), allocatable                  :: periods
  integer                                             :: rank, nprocs, err, comm3d
  integer                                             :: left, right, down, top
  integer                                             :: front, back
  !!---------------PARAMETERS FOR COOLING-----------------!!
  real, dimension(:,:), allocatable                   :: cooltable
  real                                                :: cool_Tmax, cool_Tmin
  integer                                             :: nptsT
  !!------------------ONLY FOR WINDS----------------------!!
!#ifdef PGRV
  real, dimension(:,:), allocatable                   :: posstars
  integer                                             :: Nstrs
!#endif
  
end module globals
  
