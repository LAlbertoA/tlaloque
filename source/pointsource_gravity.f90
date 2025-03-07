module pgravity

#ifdef GRAV
  
  use constants
  use globals, only: posstars, Nstrs, phip
  
  implicit none

contains
  
  subroutine pointmass_potential(phi)

    use parameters, only: nx, ny, nz

    implicit none

    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout)     :: phi
    real, dimension(0:Nstrs-1)                               :: starsmass
    integer                                                  :: i, j, k, s
    real                                                     :: x, y, z, r

    starsmass(:) = 8.0*MSUN
    
    do i = 0,nx+1
       call xcoord(i,x)
       do j = 0,ny+1
          call ycoord(j,y)
          do k = 0,nz+1
             call zcoord(k,z)
             do s = 0,Nstrs-1
                r = sqrt((posstars(0,s)-x)**2 + (posstars(1,s)-y)**2 + (posstars(2,s)-z)**2)
                phi(i,j,k) = phi(i,j,k) + GR*starsmass(s)/r
             enddo
          enddo
       enddo
    enddo

  end subroutine pointmass_potential

#endif
  
end module pgravity
