module user

  use parameters
  use globals, only: posstars, Nstrs
#ifdef GRAV
  use globals, only: PHIP
#endif
  use winds
  use pgravity

  implicit none

  type(spherical_wind_type), dimension(:), allocatable :: spherical_wind

contains

  !! Subroutine to set the initial conditions of the problem to be modeled.
  !! If using the winds module, set the parameters and make the call to the
  !! `imposeSphericalWind` subroutine of each star. The positions of stars
  !! read from the `posfile` are stored in the array `posstars` already in
  !! scope of this module
  subroutine user_initconds(U)

    use constants
    use parameters, only: dx
    
    implicit none

    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(inout) :: U
    integer                                                                    :: i

!    call load_objects()
    
!    do i = 0,Nstrs-1
       
!       spherical_wind(i)%xc = posstars(0,i)
!       spherical_wind(i)%yc = posstars(1,i)
!       spherical_wind(i)%zc = posstars(2,i)
!       spherical_wind(i)%radius = 0.01876*PC
!       spherical_wind(i)%mdot = 1.0e-6 * MSUN/YEAR
!       spherical_wind(i)%vinf = 200.0 * KMS
!       spherical_wind(i)%temp = 5000.0
!       spherical_wind(i)%mu = mu0
!       
!       call imposeSphericalWind(spherical_wind(i), U)

!    enddo
    call ExplosionPuntual()
#ifdef GRAV
    call pointmass_potential(PHIP)
#endif     
  end subroutine user_initconds

  !! Subroutine to Set the boundary conditions of the problem.
  !! If using the winds module, place here the call to `imposeSphericalWind`
  !! subroutine to each star.
  subroutine user_boundaries(U)

    implicit none

    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(inout) :: U
    integer                                                                    :: i

    do i = 0,Nstrs-1

       call imposeSphericalWind(spherical_wind(i), U)

    enddo

  end subroutine user_boundaries

  !! Subroutine to load star positions for the winds and pointsource_gravity module
  !! First line of file should be the number of stars. From second line to last line
  !! should be the positions of the stars, only one star per line and ordered like: `x  y  z`
  !! with blank space between coordinates
  subroutine load_objects()

    implicit none

    integer                           :: i

    open(unit=33,file=posfile,status="old",action="read")

    read(33,*) Nstrs
    
    allocate(spherical_wind(0:Nstrs-1))
    allocate(posstars(0:2,0:Nstrs-1))
    
    do i = 0,Nstrs-1
       read(33,*) posstars(:,i)
    enddo

    close(33)

  end subroutine load_objects

  !! Example problems 
  subroutine ExplosionPuntual()

    use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma, xr, yr, zr
    use globals
    use constants
    
    implicit none
    
    real, parameter       :: rhoO = 0.125, uO = 0.0, vO = 0.0, wO = 0.0, pO = 0.1
    real, parameter       :: rhoI = 1.0, uI = 0.0, vI = 0.0, wI = 0.0, pIn = 1.0
    real, parameter       :: xc = 0.0, yc = 0.0, zc = 0.0, rad = 0.003*PC
    integer               :: i, j, k
    real                  :: x, y, z, r
    
    x = 0.0
    y = 0.0
    z = 0.0

    do i = nxmin, nxmax
       do j = nymin, nymax
          do k = nzmin, nzmax
             call xcoord(i, x)
             call ycoord(j, y)
             call zcoord(k, z)
             r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
             if (r > 0.4) then
                U(1,i,j,k) = rhoO
                U(2,i,j,k) = rhoO*uO
                U(3,i,j,k) = rhoO*vO
                U(4,i,j,k) = rhoO*wO
                U(5,i,j,k) = 0.5*rhoO*(uO**2 + vO**2 + wO**2)+pO/(gamma-1)
             else
                U(1,i,j,k) = rhoI
                U(2,i,j,k) = rhoI*uI
                U(3,i,j,k) = rhoI*vI
                U(4,i,j,k) = rhoI*wI
                U(5,i,j,k) = 0.5*rhoI*(uI**2 + vI**2 + wI**2)+pIn/(gamma-1)
             endif
          enddo
       enddo
    enddo
#ifdef GRAV
    PHI(:,:,:) = 0.0
#endif
  end subroutine ExplosionPuntual
  
  subroutine JeansColapse(lambda)
    
    use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma, xr, yr, zr
    use globals
    use constants
    
    implicit none
    
    real, intent(in)      :: lambda
    real, parameter       :: uO = 0.0, vO = 0.0, wO = 0.0, rhoO = 1.0
    real, parameter       :: uI = 0.0, vI = 0.0, wI = 0.0
    real, parameter       :: xc = 0, yc = 0, zc = 0
    integer               :: i, j, k
    real                  :: x, y, z, pO, rho
    
    x = 0.0
    y = 0.0
    z = 0.0
    
    do i = nxmin, nxmax
       do j = nymin, nymax
          do k = nzmin, nzmax
             call xcoord(i, x)
             call ycoord(j, y)
             call zcoord(k, z)
             rho = rhoO*(1.0 + 0.01*cos(2.0*pi*x/2.0))
             pO = ((lambda**2)/pi)*rho
             U(1,i,j,k) = rho
             U(2,i,j,k) = 0.0
             U(3,i,j,k) = 0.0
             U(4,i,j,k) = 0.0
             U(5,i,j,k) = 0.5*rhoO*(uO**2 + vO**2 + wO**2)+pO/(gamma-1)
          enddo
       enddo
    enddo
#ifdef GRAV
    PHI(:,:,:) = 0.0
#endif
  end subroutine JeansColapse
  
  subroutine EvrardCollapse()
    
    use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma
    use globals
    use constants
    
    implicit none
    
    real, parameter       :: uO = 0.0, vO = 0.0, wO = 0.0, pO = (2.0/3.0)*0.05
    real, parameter       :: uI = 0.0, vI = 0.0, wI = 0.0
    integer               :: i, j, k
    real                  :: x, y, z, r, M, R1
    
    x = 0.0
    y = 0.0
    z = 0.0
    
    M = 1.0
    R1 = 1.0
    
    do i = nxmin, nxmax
       do j = nymin, nymax
          do k = nzmin, nzmax
             call xcoord(i, x)
             call ycoord(j, y)
             call zcoord(k, z)
             r = sqrt(x**2+y**2+z**2)
             if (r>R1) then
                U(1,i,j,k) = 1.e-4
                U(2,i,j,k) = 0.0
                U(3,i,j,k) = 0.0
                U(4,i,j,k) = 0.0
                U(5,i,j,k) = 1.e-4*0.05
             else
                U(1,i,j,k) = M/(2.0*PI*r*R1**2)
                U(2,i,j,k) = 0.0
                U(3,i,j,k) = 0.0
                U(4,i,j,k) = 0.0
                U(5,i,j,k) = 0.05/(2.0*PI*r*R1**2)
             end if
          enddo
       enddo
    enddo
#ifdef GRAV
    PHI(:,:,:) = 0.0
#endif
  end subroutine EvrardCollapse

  subroutine NonRotatingSphere()

    use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma, xr, yr, zr, dx
    use globals
    use constants

    implicit none

    real, parameter       :: xc = 0, yc = 0, zc = 0, cs = 0.167*KMS
    real, parameter       :: uO = 0.0, vO = 0.0, wO = 0.0, rhoO = 0.01*1.0e-15, pO = (rhoO*cs**2)
    real, parameter       :: uI = 0.0, vI = 0.0, wI = 0.0, rhoI = 1.0e-15, pIn = pO
    integer               :: i, j, k
    real                  :: x, y, z, r, R1

    x = 0.0
    y = 0.0
    z = 0.0

    R1 = 7.8e15

    do i = nxmin, nxmax
       do j = nymin, nymax
          do k = nzmin, nzmax
             call xcoord(i, x)
             call ycoord(j, y)
             call zcoord(k, z)
             r = sqrt(x**2+y**2+z**2)
             if (r>R1) then
                U(1,i,j,k) = rhoO
                U(2,i,j,k) = 0.0
                U(3,i,j,k) = 0.0
                U(4,i,j,k) = 0.0
                U(5,i,j,k) = pO/(gamma-1)
             else
                U(1,i,j,k) = rhoI
                U(2,i,j,k) = 0.0
                U(3,i,j,k) = 0.0
                U(4,i,j,k) = 0.0
                U(5,i,j,k) = pIn/(gamma-1)
             end if
          enddo
       enddo
    enddo
#ifdef GRAV
    PHI(:,:,:) = 0.0
#endif
  end subroutine NonRotatingSphere

  subroutine RCW120()

    use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma, dx, mu0
    use globals
    use constants
    use winds

    implicit none

    type(spherical_wind_type) :: spherical_wind
    type(plane_wind_type)     :: plane_wind
    real, parameter           :: pO = 3000*KB*100
    real, parameter           :: uO = 5.0*KMS, vO = 0.0, wO = 0.0, rhoO = 3000*mu0*AMU
    real, parameter           :: uI = 5.0*KMS, vI = 0.0, wI = 0.0, rhoI = 1.0e-15
    real                      :: x, y, z, r
    integer                   :: i, j, k

    plane_wind%plane = 1
    plane_wind%rho = rhoO
    plane_wind%vel = 5.0*KMS
    plane_wind%temp = 100.0
    plane_wind%mu = mu0

    spherical_wind%xc = 0.5*PC
    spherical_wind%yc = 1.5*PC
    spherical_wind%zc = 1.5*PC
    spherical_wind%radius = 0.09*PC
    spherical_wind%mdot = 2.7e-7 * MSUN/YEAR
    spherical_wind%vinf = 1000.0 * KMS
    spherical_wind%temp = 1.0e6
    spherical_wind%mu = mu0

    do i = nxmin, nxmax
       do j = nymin, nymax
          do k = nzmin, nzmax
             U(1,i,j,k) = rhoO
             U(2,i,j,k) = rhoO*uO
             U(3,i,j,k) = rhoO*vO
             U(4,i,j,k) = rhoO*wO
             U(5,i,j,k) = 0.5*rhoO*(uO**2 + vO**2 + wO**2) + pO/(gamma-1)
          enddo
       enddo
    enddo
    call imposeSphericalWind(spherical_wind, U)
    call imposePlaneWind(plane_wind, U)
#ifdef GRAV
    PHI(:,:,:) = 0.0
    do i = 0, nxmax-1
       do j = 0, nymax-1
          call xcoord(i,x)
          call ycoord(j,y)
          k = 0
          call zcoord(k,z)
          r = sqrt(x**2+y**2+z**2)
          PHI(i,j,k) = -GR*(rhoO*(4.0/3.0)*PI*r**3)/r
          k = nzmax-1
          call zcoord(k,z)
          r = sqrt(x**2+y**2+z**2)
          PHI(i,j,k) = -GR*(rhoO*(4.0/3.0)*PI*r**3)/r
       enddo
    enddo

    do j = 0, nymax-1
       do k = 0, nzmax-1
          call ycoord(j,y)
          call zcoord(k,z)
          i = 0
          call xcoord(i,x)
          r = sqrt(x**2+y**2+z**2)
          PHI(i,j,k) = -GR*(rhoO*(4.0/3.0)*PI*r**3)/r
          i = nxmax-1
          call xcoord(i,x)
          r = sqrt(x**2+y**2+z**2)
          PHI(i,j,k) = -GR*(rhoO*(4.0/3.0)*PI*r**3)/r
       enddo
    enddo

    do i = 0, nxmax-1
       do k = 0, nzmax-1
          call xcoord(i,x)
          call zcoord(k,z)
          j = 0
          call ycoord(j,y)
          r = sqrt(x**2+y**2+z**2)
          PHI(i,j,k) = -GR*(rhoO*(4.0/3.0)*PI*r**3)/r
          j = nymax-1
          call ycoord(j,y)
          r = sqrt(x**2+y**2+z**2)
          PHI(i,j,k) = -GR*(rhoO*(4.0/3.0)*PI*r**3)/r
       enddo
    enddo
#endif
    
  end subroutine RCW120

  subroutine Stars()

    use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma, mu0
    use globals
    use constants
    use winds

    implicit none

    real, parameter           :: pO = 1.0e4*mu0*AMU*KB*100/(mu0*AMU)
    real, parameter           :: uO = 0.0, vO = 0.0, wO = 0.0, rhoO = 1.0e4*mu0*AMU
    real, parameter           :: uI = 0.0, vI = 0.0, wI = 0.0, rhoI = 1.0e-15
    integer                   :: i, j, k

    do i = nxmin, nxmax
       do j = nymin, nymax
          do k = nzmin, nzmax
             U(1,i,j,k) = rhoO
             U(2,i,j,k) = rhoO*uO
             U(3,i,j,k) = rhoO*vO
             U(4,i,j,k) = rhoO*wO
             U(5,i,j,k) = 0.5*rhoO*(uO**2 + vO**2 + wO**2) + pO/(gamma-1)
          enddo
       enddo
    enddo

#ifdef GRAV
    PHI(:,:,:) = 0.0
#endif

  end subroutine Stars

end module user
