module constants

  implicit none

  ! Fundamental and astrophysical constants (cgs)
  real, parameter :: AMU  = 1.660538782e-24   ! Atomic Mass Unit
  real, parameter :: KB   = 1.380650400e-16   ! Boltzmann constant
  real, parameter :: GR   = 6.67430e-8        ! Gravitational constant
  real, parameter :: PC   = 3.085677588e+18   ! Parsec
  real, parameter :: AU   = 1.495978707e+13   ! Astronomical unit
  real, parameter :: YEAR = 3.155673600e+7    ! Year (Earth, sidereal)
  real, parameter :: KYR  = 3.155673600e+10   ! One thousand years
  real, parameter :: MSUN = 1.988920000e+33   ! Solar mass
  real, parameter :: KMS  = 1.0e5             ! km/s in cgs
  real, parameter :: PI   = 3.14159265358979  ! Ratio of perimeter to diameter

  ! Named Constants

  ! BC types
  integer, parameter :: BC_OUTFLOW    = 1      ! Free flow (zero gradient) BC
  integer, parameter :: BC_REFLECTIVE = 2      ! Reflection (velocity) BC
  integer, parameter :: BC_PERIODIC   = 3      ! Periodic BC
  integer, parameter :: BC_INFLOW     = 4      ! Inflow BCs

  ! Slope limiters
  integer, parameter :: LIMITER_NONE     = 1   ! No limiter
  integer, parameter :: LIMITER_VANLEER  = 2   ! Van Leer limiter
  integer, parameter :: LIMITER_MINMOD   = 3   ! Minmod limiter
  integer, parameter :: LIMITER_UMIST    = 4   ! UMIST limiter
  integer, parameter :: LIMITER_WOODWARD = 5   ! Woodoward limiter
  integer, parameter :: LIMITER_SUPERBEE = 6   ! Superbee limiter

  ! Output file type

  integer, parameter :: DAT = 1
  integer, parameter :: VTK = 2
  integer, parameter :: BIN = 3
  
end module
