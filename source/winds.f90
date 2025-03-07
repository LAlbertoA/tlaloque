module winds

  implicit none

  type spherical_wind_type
     real :: xc, yc, zc
     real :: vx, vy, vz
     real :: radius
     real :: mdot
     real :: vinf
     real :: temp
     real :: mu
  end type spherical_wind_type

  integer, parameter           :: PLANE_LEFT   = 1
  integer, parameter           :: PLANE_RIGHT  = 2
  integer, parameter           :: PLANE_FRONT  = 3
  integer, parameter           :: PLANE_BACK   = 4
  integer, parameter           :: PLANE_BOTTOM = 5
  integer, parameter           :: PLANE_TOP    = 6
  
  type plane_wind_type
     integer :: plane
     real :: rho
     real :: vel
     real :: temp
     real :: mu
  end type plane_wind_type
  
contains

  subroutine imposeSphericalWind (wind_params, uvars)

    use parameters
    use globals

    implicit none

    type(spherical_wind_type), intent(in) :: wind_params
    real, intent(inout)                   :: uvars(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer                               :: i, j, k
    real                                  :: xc, yc, zc, vwx, vwy, vwz, radius, mdot, vinf, temp
    real                                  :: dens, vx, vy, vz, pres, x, y, z, dist, mu

    ! Unpack wind source parameters
    xc = wind_params%xc
    yc = wind_params%yc
    zc = wind_params%zc
    vwx = wind_params%vx
    vwy = wind_params%vy
    vwz = wind_params%vz
    radius = wind_params%radius
    mdot = wind_params%mdot
    vinf = wind_params%vinf
    temp = wind_params%temp
    mu = wind_params%mu

    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nzmin,nzmax
             
             call xcoord(i,x)
             call ycoord(j,y)
             call zcoord(k,z)
             
             dist = sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)   ! phys units
             
             if (dist.lt.radius) then

                ! Calculate wind density, velocity and pressure in this cell
                dens = mdot/(vinf*dist*dist*4.0*PI)
                vx = vinf*(x-xc)/dist + vwx
                vy = vinf*(y-yc)/dist + vwy
                vz = vinf*(z-zc)/dist + vwz
                pres= dens*KB*temp/(mu*AMU)
                
                ! Convert primitives and set flow vars for this cell
                uvars(1,i,j,k) = dens
                uvars(2,i,j,k) = dens*vx
                uvars(3,i,j,k) = dens*vy
                uvars(4,i,j,k) = dens*vz
                uvars(5,i,j,k) = 0.5*dens*(vx**2 + vy**2 + vz**2) + pres/(gamma-1)
             end if
             
          end do
       end do
    end do
    
  end subroutine imposeSphericalWind

  subroutine imposePlaneWind(wind_params, uvars)

    use parameters
    use globals
    implicit none

    type(plane_wind_type), intent(in) :: wind_params
    real, intent(inout) :: uvars(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: i, j, k, plane
    real :: dens, vx, vy, vz, pres, temp, vel, mu

    ! Unpack wind source parameters
    plane = wind_params%plane
    dens  = wind_params%rho
    vel   = wind_params%vel
    temp  = wind_params%temp
    mu    = wind_params%mu

    do i=nxmin,nxmax
       do j=nymin,nymax
          do k=nz,nzmax
             
             ! Compute velocity components and pressure
             vx = 0.0
             vy = 0.0
             vz = 0.0
             select case(plane)
             case (PLANE_LEFT)
                vx = vel
             case (PLANE_RIGHT)
                vx = -vel
             case (PLANE_FRONT)
                vy = vel
             case (PLANE_BACK)
                vy = -vel
             case (PLANE_BOTTOM)
                vz = vel
             case (PLANE_TOP)
                vz = -vel
             end select
             pres = dens/(mu*AMU)*KB*temp
             
             ! Convert primitives and set flow vars for this cell
             
             uvars(1,i,j,k) = dens
             uvars(2,i,j,k) = dens*vx
             uvars(3,i,j,k) = dens*vy
             uvars(4,i,j,k) = dens*vz
             uvars(5,i,j,k) = 0.5*dens*(vx**2 + vy**2 + vz**2) + pres/(gamma-1)
             
          end do
       end do
    end do

  end subroutine imposePlaneWind

end module winds
