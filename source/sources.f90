module sources

  implicit none

contains
  
#ifdef GRAV

  real function Legendre(x,l)

    implicit none

    real, intent(in)    :: x
    integer, intent(in) :: l

    select case(l)

    case(0)
       Legendre = 1.0
    case(1)
       Legendre = x
    case(2)
       Legendre = 0.5*(3.0*x**2 - 1.0)
    case(3)
       Legendre = 0.5*(5.0*x**3 - 3.0*x)
    case(4)
       Legendre = (1.0/8.0)*(35.0*x**4 - 30.0*x**2 + 3.0)
    end select

  end function Legendre

  real function aLegendre(x,l,m)

    implicit none

    real, intent(in)     :: x
    integer, intent(in)  :: l, m

    select case(l)

    case(0)
       aLegendre = 1.0
    case(1)

       select case(m)

       case(0)
          aLegendre = x
       case(1)
          aLegendre = -(1.0 - x**2)**(1.0/2.0)
       end select

    case(2)

       select case(m)

       case(0)
          aLegendre = 0.50*(3.0*x**2 - 1.0)
       case(1)
          aLegendre = -3.0*x*(1.0 - x**2)**(1.0/2.0)
       case(2)
          aLegendre = 3.0*(1.0 - x**2)
       end select

    case(3)

       select case(m)

       case(0)
          aLegendre = 0.50*(5.0*x**3 - 3.0*x)
       case(1)
          aLegendre = (3.0/2.0)*(1.0 - 5.0*x**2)*(1.0 - x**2)**(1.0/2.0)
       case(2)
          aLegendre = 15.0*x*(1.0 - x**2)
       case(3)
          aLegendre = -15.0*(1.0 - x**2)**(3.0/2.0)
       end select

    case(4)

       select case(m)

       case(0)
          aLegendre = (1.0/8.0)*(35.0*x**4 - 30.0*x**2 + 3.0)
       case(1)
          aLegendre = -(5.0/2.0)*(7.0*x**3 - 3.0*x)*(1.0 - x**2)**(1.0/2.0)
       case(2)
          aLegendre = (15.0/2.0)*(7.0*x**2 - 1.0)*(1.0 - x**2)
       case(3)
          aLegendre = -105.0*x*(1.0 - x**2)**(3.0/2.0)
       case(4)
          aLegendre = 105.0*(1.0 - x**2)**2
       end select

    end select

  end function aLegendre
  
  subroutine MultiGrid(dens, phi)

    use globals, only: Error, Residue, ret, rank, left, right, down, top, front, back, it
    use parameters, only: nx, ny, nz, lvlm, logged, logu, Gconst
    use constants
    
    implicit none

    real, intent(in)                     :: dens(0:nx+1,0:ny+1,0:nz+1)
    real, intent(inout)                  :: phi(0:nx+1,0:ny+1,0:nz+1)
    real, dimension(4,4)                 :: cosMoment, sinMoment
    real, dimension(0:4)                 :: c0
    integer, parameter                   :: Solve = 1, Relax = 2
    integer                              :: nxl, nyl, nzl, lv, level, iteraciones, iter
    real                                 :: Dif

    if (mod(it,10).eq.0) then
       call Moments(c0,cosMoment,sinMoment,dens)
       if (left<0) then
          call phi_xboundaries(0,c0,cosMoment,sinMoment,phi)
       endif
       if (right<0) then
          call phi_xboundaries(nx+1,c0,cosMoment,sinMoment,phi)
       endif
       if (down<0) then
          call phi_yboundaries(0,c0,cosMoment,sinMoment,phi)
       endif
       if (top<0) then
          call phi_yboundaries(ny+1,c0,cosMoment,sinMoment,phi)
       endif
       if (front<0) then
          call phi_zboundaries(0,c0,cosMoment,sinMoment,phi)
       endif
       if (back<0) then
          call phi_zboundaries(nz+1,c0,cosMoment,sinMoment,phi)
       endif
    endif
    Residue(0)%data = 4.0*PI*Gconst*dens
    Error(0)%data = phi
    Dif = 1.e28
    iter = 0
    do while (Dif > 10.0 .and. iter < 31)
       iter = iter + 1
!    do iter = 1,10
       iteraciones = 6

!       print*, "Finest Grid ------------------------------------------"

       do lv = 1,lvlm
          Error(lv)%data = 0.0
       enddo

       level = 0
       
       call poisson_solver(Residue(0)%data,Error(0)%data,nx,ny,nz,Relax,iteraciones,Dif,level)

       do level = 1,lvlm

          nxl = int(nx/2**level)
          nyl = int(ny/2**level)
          nzl = int(nz/2**level)

!          print*, "Restriction on ", level-1

          if (level == lvlm) then
             call resid(nxl*2,nyl*2,nzl*2,Error(level-1)%data,Residue(level-1)%data,ret(level-1)%data)
             call restriction(nxl*2,nyl*2,nzl*2,nxl,nyl,nzl,ret(level-1)%data,Residue(level)%data)
!             print*, "Solving residue equation on level", level
             call poisson_solver(Residue(level)%data,Error(level)%data,nxl,nyl,nzl,Solve,iteraciones,Dif,level)
          else
             call resid(nxl*2,nyl*2,nzl*2,Error(level-1)%data,Residue(level-1)%data,ret(level-1)%data)
             call restriction(nxl*2,nyl*2,nzl*2,nxl,nyl,nzl,ret(level-1)%data,Residue(level)%data)
!             print*, "Relaxing residue equation on level", level
             call poisson_solver(Residue(level)%data,Error(level)%data,nxl,nyl,nzl,Relax,iteraciones,Dif,level)

          endif

       enddo

       do level = lvlm-1, 0, -1

          nxl = int(nx/2**level)
          nyl = int(ny/2**level)
          nzl = int(nz/2**level)

          ret(level)%data = 0.0

          call prolongation(int(nxl/2.0),int(nyl/2.0),int(nzl/2.0),nxl,nyl,nzl,Error(level+1)%data,ret(level)%data)

!          print*, "Prolongation on ", level

          Error(level)%data = Error(level)%data + ret(level)%data

!          print*, "Relaxing residue equation on level ", level-1
          call poisson_solver(Residue(level)%data,Error(level)%data,nxl,nyl,nzl,Relax,iteraciones,Dif,level)

       enddo

       if (rank == 0 .and. logged == .true.) then     
          write(logu,*) "Dif = ", Dif, "V-Cycle Iteration = ", iter
       endif

    enddo
    
    phi = Error(0)%data
       
  end subroutine MultiGrid
  

  subroutine  restriction(nx,ny,nz,nxl,nyl,nzl,var,varTemp)

    implicit none

    integer, intent(in)                   :: nx, ny, nz, nxl, nyl, nzl
    real, intent(inout)                   :: var(0:nx+1,0:ny+1,0:nz+1), varTemp(0:nxl+1,0:nyl+1,0:nzl+1)
    integer                               :: i, j, k, i1, j1, k1

    do i = 1,nxl
       do j = 1,nyl
          do k = 1,nzl
             i1 = i*2-1
             j1 = j*2-1
             k1 = k*2-1
             varTemp(i,j,k) = (1.0/8.0)*(var(i1,j1,k1)+var(i1+1,j1,k1)+var(i1,j1+1,k1)+var(i1,j1,k1+1)+var(i1+1,j1+1,k1)+&
                  var(i1+1,j1,k1+1)+var(i1,j1+1,k1+1)+var(i1+1,j1+1,k1+1))
          end do
       end do
    end do

  end subroutine restriction
  
  subroutine prolongation(nxl,nyl,nzl,nx,ny,nz,var,varTemp)

    implicit none

    integer, intent(in)                   :: nx, ny, nz, nxl, nyl, nzl
    real, intent(inout)                   :: var(0:nxl+1,0:nyl+1,0:nzl+1), varTemp(0:nx+1,0:ny+1,0:nz+1)
    integer                               :: i, j, k, i1, j1, k1

    do i = 1,nxl
       do j = 1,nyl
          do k = 1,nzl
             i1 = i*2-1
             j1 = j*2-1
             k1 = k*2-1
             varTemp(i1,j1,k1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i-1,j,k)+9*var(i,j-1,k)+9*var(i,j,k-1)+&
                  3*var(i-1,j-1,k)+3*var(i-1,j,k-1)+3*var(i,j-1,k-1)+var(i-1,j-1,k-1))

             varTemp(i1+1,j1,k1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i+1,j,k)+9*var(i,j-1,k)+9*var(i,j,k-1)+&
                  3*var(i+1,j-1,k)+3*var(i+1,j,k-1)+3*var(i,j-1,k-1)+var(i+1,j-1,k-1))

             varTemp(i1,j1+1,k1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i-1,j,k)+9*var(i,j+1,k)+9*var(i,j,k-1)+&
                  3*var(i-1,j+1,k)+3*var(i-1,j,k-1)+3*var(i,j+1,k-1)+var(i-1,j+1,k-1))

             varTemp(i1+1,j1+1,k1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i+1,j,k)+9*var(i,j+1,k)+9*var(i,j,k-1)+&
                  3*var(i+1,j+1,k)+3*var(i+1,j,k-1)+3*var(i,j+1,k-1)+var(i+1,j+1,k-1))

             varTemp(i1,j1,k1+1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i-1,j,k)+9*var(i,j-1,k)+9*var(i,j,k+1)+&
                  3*var(i-1,j-1,k)+3*var(i-1,j,k+1)+3*var(i,j-1,k+1)+var(i-1,j-1,k+1))

             varTemp(i1+1,j1,k1+1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i+1,j,k)+9*var(i,j-1,k)+9*var(i,j,k+1)+&
                  3*var(i+1,j-1,k)+3*var(i+1,j,k+1)+3*var(i,j-1,k+1)+var(i+1,j-1,k+1))

             varTemp(i1,j1+1,k1+1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i-1,j,k)+9*var(i,j+1,k)+9*var(i,j,k+1)+&
                  3*var(i-1,j+1,k)+3*var(i-1,j,k+1)+3*var(i,j+1,k+1)+var(i-1,j+1,k+1))

             varTemp(i1+1,j1+1,k1+1) = (1.0/64.0)*(27*var(i,j,k)+9*var(i+1,j,k)+9*var(i,j+1,k)+9*var(i,j,k+1)+&
                  3*var(i+1,j+1,k)+3*var(i+1,j,k+1)+3*var(i,j+1,k+1)+var(i+1,j+1,k+1))
          enddo
       enddo
    enddo

  end subroutine prolongation

  subroutine resid(nx,ny,nz,var,source,resi)

    use parameters, only: xl, xr, mpix
    implicit none

    integer, intent(in)   :: nx, ny, nz
    real, intent(in)      :: var(0:nx+1,0:ny+1,0:nz+1), source(0:nx+1,0:ny+1,0:nz+1)
    real, intent(out)     :: resi(0:nx+1,0:ny+1,0:nz+1)
    real                  :: grad(0:nx+1,0:ny+1,0:nz+1)
    real                  :: dxl
    integer               :: i, j, k

    grad = 0.0

    dxl = (xr-xl)/(mpix*nx)

    do i = 1, nx
       do j = 1, ny
          do k = 1, nz
             grad(i,j,k) = (var(i-1,j,k)+var(i+1,j,k)+var(i,j-1,k)+var(i,j+1,k)+var(i,j,k-1)+var(i,j,k+1)-6.0*var(i,j,k))/(dxl**2)
          enddo
       enddo
    enddo

    !  call laplacian(nx,ny,nz,var,grad)
    resi = source - grad

  end subroutine resid

  subroutine  poisson_solver (source,var,nx,ny,nz,relsolv,it,Dif,level)

    use parameters, only: xl, xr, size, mpi_real_kind, mpix
#ifdef MPIP
    use globals, only: rank, nprocs, err, left, right, down, top, front, back, comm3d
#endif
    implicit  none

#ifdef MPIP
    include "mpif.h"
#endif
    
    integer, intent(in)       ::  nx, ny, nz, relsolv, it, level
    real, intent(in)          ::  source(0:nx+1,0:ny+1,0:nz+1)
    real, intent(inout)       ::  var(0:nx+1,0:ny+1,0:nz+1)
    real, intent(inout)       ::  Dif
    real, parameter           ::  gam = 1.      !! gam = 1 -> Gauss-Seidel, 1 < gam < 2 -> SOR
    integer, parameter        ::  Solve = 1, Relax = 2

    real                      ::  dxl, er, d1, difp, w
    integer                   ::  i, j, k, iter, start

    dxl = (xr-xl)/(mpix*nx)
    d1 = 0.0
    
    select case (relsolv)
       
    case(Solve)
       
       w = 1.0
       Dif = 1000.0
       iter = 0
       er = 0.001

       if (level>0) then
#ifdef MPIP
          if (left<0) then
             var(0,:,:) = -var(1,:,:)
          endif
          if (right<0) then
             var(nx+1,:,:) = -var(nx,:,:)
          endif
          if (down<0) then
             var(:,0,:) = -var(:,1,:)
          endif
          if (top<0) then
             var(:,ny+1,:) = -var(:,ny,:)
          endif
          if (front<0) then
             var(:,:,0) = -var(:,:,1)
          endif
          if (back<0) then
             var(:,:,nz+1) = -var(:,:,nz)
          endif
#else
          var(0,:,:) = -var(1,:,:)
          var(:,0,:) = -var(:,1,:)
          var(:,:,0) = -var(:,:,1)
          var(nx+1,:,:) = -var(nx,:,:)
          var(:,ny+1,:) = -var(:,ny,:)
          var(:,:,nz+1) = -var(:,:,nz)
#endif
       endif       
       
       do while (Dif>er)

          !! Successive Over-Relaxation
          Dif = 0
          difp = 0
#ifdef MPIP
          call PotentialBoundaries(nx, ny, nz, var)
#endif
          do i = 1,nx
             do j = 1,ny
                do k = 1,nz
                   d1 = var(i,j,k)
                   var(i,j,k) = ((1.0/6.0)*(var(i+1,j,k) + var(i-1,j,k) + var(i,j+1,k) + var(i,j-1,k) + var(i,j,k+1) + &
                        var(i,j,k-1) - source(i,j,k)*dxl**2))*w + (1.0-w)*var(i,j,k)
                   difp = difp + abs(var(i,j,k)-d1)
                enddo
             enddo
          enddo

          iter = iter + 1

#ifdef MPIP
          call mpi_allreduce(difp, Dif, 1, mpi_real_kind, mpi_sum, comm3d, err)
#else
          Dif = difp
#endif
!          if (rank == 0) then
!             write (*,*) 'Dif = ', Dif, 'iter = ', iter
!          endif

       end do
       
    case(Relax)

       w = 1.0
       if (level>0) then
#ifdef MPIP
          if (left<0) then
             var(0,:,:) = -var(1,:,:)
          endif
          if (right<0) then
             var(nx+1,:,:) = -var(nx,:,:)
          endif
          if (down<0) then
             var(:,0,:) = -var(:,1,:)
          endif
          if (top<0) then
             var(:,ny+1,:) = -var(:,ny,:)
          endif
          if (front<0) then
             var(:,:,0) = -var(:,:,1)
          endif
          if (back<0) then
             var(:,:,nz+1) = -var(:,:,nz)
          endif
#else
          var(0,:,:) = -var(1,:,:)
          var(:,0,:) = -var(:,1,:)
          var(:,:,0) = -var(:,:,1)
          var(nx+1,:,:) = -var(nx,:,:)
          var(:,ny+1,:) = -var(:,ny,:)
          var(:,:,nz+1) = -var(:,:,nz)
#endif
       endif

       do iter = 1,it-1
          !! Successive Over-Relaxation
          Dif = 0
          difp = 0
#ifdef MPIP
          call PotentialBoundaries(nx, ny, nz, var)
#endif
          do i = 1,nx
             do j = 1,ny
                if (((mod(i,2).eq.0).and.(mod(j,2).eq.1)).or.((mod(i,2).eq.0).and.(mod(j,2).eq.1))) then
                   start = 1
                else
                   start = 2
                endif
                do k = start,nz,2
                   d1 = var(i,j,k)
                   var(i,j,k) = ((1.0/6.0)*(var(i+1,j,k) + var(i-1,j,k) + var(i,j+1,k) + var(i,j-1,k) + var(i,j,k+1) + &
                        var(i,j,k-1) - source(i,j,k)*dxl**2))*w + (1.0-w)*var(i,j,k)
                   difp = difp + abs(var(i,j,k)-d1)
                enddo
             enddo
          enddo
#ifdef MPIP
          call PotentialBoundaries(nx, ny, nz, var)
#endif
          do i = 1,nx
             do j = 1,ny
                if (((mod(i,2).eq.0).and.(mod(j,2).eq.1)).or.((mod(i,2).eq.0).and.(mod(j,2).eq.1))) then
                   start = 2
                else
                   start = 1
                endif
                do k = start,nz,2
                   d1 = var(i,j,k)
                   var(i,j,k) = ((1.0/6.0)*(var(i+1,j,k) + var(i-1,j,k) + var(i,j+1,k) + var(i,j-1,k) + var(i,j,k+1) + &
                        var(i,j,k-1) - source(i,j,k)*dxl**2))*w + (1.0-w)*var(i,j,k)
                   difp = difp + abs(var(i,j,k)-d1)
                enddo
             enddo
          enddo

#ifdef MPIP
          call mpi_allreduce(difp, Dif, 1, mpi_real_kind, mpi_sum, comm3d, err)
#else
          Dif = difp
#endif          
       enddo

    end select

!    if (rank == 0) then
!       write (*,*) 'Dif = ', Dif, 'iter = ', iter
!    endif
    
  end subroutine poisson_solver

#ifdef MPIP

  subroutine PotentialBoundaries(nx, ny, nz, var)

    use parameters, only: mpi_real_kind
    use globals, only: rank, nprocs, err, left, right, down, top, front, back, comm3d
    implicit none

    include "mpif.h"

    integer, intent(in)                                     :: nx, ny, nz
    real(8), intent(inout), dimension(0:nx+1,0:ny+1,0:nz+1) :: var
    real(8), dimension(0:ny+1,0:nz+1)                       :: sendr, sendl, recvr, recvl
    real(8), dimension(0:nx+1,0:nz+1)                       :: sendt, sendd, recvt, recvd
    real(8), dimension(0:nx+1,0:ny+1)                       :: sendb, sendf, recvb, recvf
    integer                                                 :: status(MPI_STATUS_SIZE)
    integer                                                 :: srsizex, srsizey, srsizez

    srsizex = (ny+2)*(nz+2)
    srsizey = (nx+2)*(nz+2)
    srsizez = (nx+2)*(ny+2)
    
    sendr(0:ny+1,0:nz+1) = var(nx,0:ny+1,0:nz+1)
    sendl(0:ny+1,0:nz+1) = var(1 ,0:ny+1,0:nz+1)
    sendt(0:nx+1,0:nz+1) = var(0:nx+1,ny,0:nz+1)
    sendd(0:nx+1,0:nz+1) = var(0:ny+1,1 ,0:nz+1)
    sendb(0:nx+1,0:ny+1) = var(0:nx+1,0:ny+1,nz)
    sendf(0:nx+1,0:ny+1) = var(0:nx+1,0:ny+1,1 )
    
    call mpi_sendrecv(sendr, srsizex, mpi_real_kind, right, 0,        &
         recvl, srsizex, mpi_real_kind, left , 0, comm3d, status, err)
    !    print'(A,i3,A,i3,A,i3)', 'Arrays sended to right:', right, 'and recieved from left:', left, ' in rank ', rank
    call mpi_sendrecv(sendl, srsizex, mpi_real_kind, left , 0,        &
         recvr, srsizex, mpi_real_kind, right, 0, comm3d, status, err)
    !    print'(A,i3,A,i3,A,i3)', 'Arrays sended to left:', left, 'and recieved from right:', right, ' in rank ', rank
    call mpi_sendrecv(sendt, srsizey, mpi_real_kind, top  , 0,        &
         recvd, srsizey, mpi_real_kind, down , 0, comm3d, status, err)
    !    print'(A,i3,A,i3,A,i3)', 'Arrays sended to top:', top, 'and recieved from down:', down, ' in rank ', rank
    call mpi_sendrecv(sendd, srsizey, mpi_real_kind, down , 0,        &
         recvt, srsizey, mpi_real_kind, top  , 0, comm3d, status, err)
    !    print'(A,i3,A,i3,A,i3)', 'Arrays sended to down:', down, 'and recieved from top:', top, ' in rank ', rank
    call mpi_sendrecv(sendb, srsizez, mpi_real_kind, back , 0,        &
         recvf, srsizez, mpi_real_kind, front, 0, comm3d, status, err)
    !    print'(A,i3,A,i3,A,i3)', 'Arrays sended to back:', back, 'and recieved from front:', front, ' in rank ', rank
    call mpi_sendrecv(sendf, srsizez, mpi_real_kind, front, 0,        &
         recvb, srsizez, mpi_real_kind, back , 0, comm3d, status, err)
    !    print'(A,i3,A,i3,A,i3)', 'Arrays sended to front:', front, 'and recieved from back:', back, ' in rank ', rank 

    if (left>=0) then
       var(0,0:ny+1,0:nz+1) = recvl(0:ny+1,0:nz+1)
    endif
    if (right>=0) then
       var(nx+1,0:ny+1,0:nz+1) = recvr(0:ny+1,0:nz+1)
    endif
    if (down>=0) then
       var(0:nx+1,0,0:nz+1) = recvd(0:nx+1,0:nz+1)
    endif
    if (top>=0) then
       var(0:nx+1,ny+1,0:nz+1) = recvt(0:nx+1,0:nz+1)
    endif
    if (front>=0) then
       var(0:nx+1,0:ny+1,0) = recvf(0:nx+1,0:ny+1)
    endif
    if (back>=0) then
       var(0:nx+1,0:ny+1,nz+1) = recvb(0:nx+1,0:ny+1)
    endif

  end subroutine PotentialBoundaries
  
#endif

  subroutine distance(r,x,y,z)

    implicit none

    real, intent(in)       :: x, y, z
    real, intent(out)      :: r

    r = sqrt(x**2 + y**2 + z**2)

  end subroutine distance

  subroutine pangle(the,x,y,z)

    implicit none

    real, intent(in)       :: x, y, z
    real, intent(out)      :: the

    the = acos(z/sqrt(x**2 + y**2 + z**2))

  end subroutine pangle

  subroutine aangle(fi,x,y)

    use constants, only: PI

    implicit none

    real, intent(in)       :: x, y
    real, intent(out)      :: fi

    fi = atan2(y,x)

    if (fi < 0) then

       fi = fi + 2.*PI

    endif

  end subroutine aangle

  subroutine Moments(c0,cm,sm,rho)

    use parameters, only: nx, ny, nz, dx, dy, dz, mpi_real_kind
    use globals, only: comm3d, err

    implicit none

    include "mpif.h"

    real, intent(in), dimension(0:nx+1,0:ny+1,0:nz+1) :: rho
    real, intent(out), dimension(4,4)                 :: cm, sm
    real, intent(out), dimension(0:4)                 :: c0
    real                                              :: r, coth, x, y, z, smoment
    real                                              :: c, cmp, smp, thet, fi, c0mp, cmoment, c0moment
    integer                                           :: i, j, k, l, m

    c0 = 0.0
    cm = 0.0
    sm = 0.0

    do l = 1,4
       do m = 1,l
          c0mp = 0.0
          cmp = 0.0
          smp = 0.0
          c = 2.0*Gamma(l-m+1.0)/Gamma(l+m+1.0)
          do i = 1,nx
             do j = 1,ny
                do k = 1,nz
                   call xcoord(i,x)
                   call ycoord(j,y)
                   call zcoord(k,z)
                   call pangle(thet,x,y,z)
                   call aangle(fi,x,y)
                   call distance(r,x,y,z)
                   coth = cos(thet)
                   c0mp = c0mp + Legendre(coth,l)*(r**l)*rho(i,j,k)*dx*dy*dz
                   cmp = cmp + aLegendre(coth,l,m)*cos(m*fi)*(r**l)*rho(i,j,k)*dx*dy*dz
                   smp = smp + aLegendre(coth,l,m)*sin(m*fi)*(r**l)*rho(i,j,k)*dx*dy*dz
                   !             print*, cmoment, smoment
                enddo
             enddo
          enddo
          call mpi_allreduce(cmp, cmoment, 1, mpi_real_kind, mpi_sum, comm3d, err)
          call mpi_allreduce(smp, smoment, 1, mpi_real_kind, mpi_sum, comm3d, err)
          cm(l,m) = c*cmoment
          sm(l,m) = c*smoment
       enddo
       call mpi_allreduce(c0mp, c0moment, 1, mpi_real_kind, mpi_sum, comm3d, err)
       c0(l) = c0moment
    enddo
    c0mp = 0.0
    l = 0
    m = 0
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             call xcoord(i,x)
             call ycoord(j,y)
             call zcoord(k,z)
             call pangle(thet,x,y,z)
             call aangle(fi,x,y)
             call distance(r,x,y,z)
             coth = cos(thet)
             c0mp = c0mp + Legendre(coth,l)*(r**l)*rho(i,j,k)*dx*dy*dz
          enddo
       enddo
    enddo
    call mpi_allreduce(c0mp, c0moment, 1, mpi_real_kind, mpi_sum, comm3d, err)
    c0(l) = c0moment

  end subroutine Moments

  subroutine phi_xboundaries(ibound,c0,cosMoment,sinMoment,phiarr)

    use parameters, only: nx, ny, nz, dx, dy, dz, Gconst

    implicit none

    integer, intent(in)                                :: ibound
    real, intent(in), dimension(0:4)                   :: c0
    real, intent(in), dimension(4,4)                   :: cosMoment, sinMoment
    real, intent(out), dimension(0:nx+1,0:ny+1,0:nz+1) :: phiarr
    integer                                            :: i, j, k, l, m
    real                                               :: phi, x, y, z, thet, fi, rb, costheta

    do j = 1,ny
       do k = 1,nz
          i = ibound
          phi = 0.0
          call xcoord(i,x)
          call ycoord(j,y)
          call zcoord(k,z)
          call pangle(thet,x,y,z)
          call aangle(fi,x,y)
          call distance(rb,x,y,z)
          costheta = cos(thet)
          phi = phi - c0(0)/(rb)
          do l = 1,4
             phi = phi - c0(l)*Legendre(costheta,l)/(rb**(l+1))
             do m = 1,l
                phi = phi - 2.0*(cosMoment(l,m)*cos(m*fi) + sinMoment(l,m)*sin(m*fi))*aLegendre(costheta,l,m)/(rb**(l+1))
             enddo
          enddo
          phiarr(i,j,k) = Gconst*phi
       enddo
    enddo

  end subroutine phi_xboundaries

  subroutine phi_yboundaries(jbound,c0,cosMoment,sinMoment,phiarr)

    use parameters, only: nx, ny, nz, dx, dy, dz, Gconst
    
    implicit none

    integer, intent(in)                                :: jbound
    real, intent(in), dimension(0:4)                   :: c0
    real, intent(in), dimension(4,4)                   :: cosMoment, sinMoment
    real, intent(out), dimension(0:nx+1,0:ny+1,0:nz+1) :: phiarr
    integer                                            :: i, j, k, l, m
    real                                               :: phi, x, y, z, thet, fi, rb, costheta

    do i = 1,nx
       do k = 1,nz
          j = jbound
          phi = 0.0
          call xcoord(i,x)
          call ycoord(j,y)
          call zcoord(k,z)
          call pangle(thet,x,y,z)
          call aangle(fi,x,y)
          call distance(rb,x,y,z)
          costheta = cos(thet)
          phi = phi - c0(0)/(rb)
          do l = 1,4
             phi = phi - c0(l)*Legendre(costheta,l)/(rb**(l+1))
             do m = 1,l
                phi = phi - 2.0*(cosMoment(l,m)*cos(m*fi) + sinMoment(l,m)*sin(m*fi))*aLegendre(costheta,l,m)/(rb**(l+1))
             enddo
          enddo
          phiarr(i,j,k) = Gconst*phi
       enddo
    enddo

  end subroutine phi_yboundaries

  subroutine phi_zboundaries(kbound,c0,cosMoment,sinMoment,phiarr)

    use parameters, only: nx, ny, nz, dx, dy, dz, Gconst

    implicit none

    integer, intent(in)                                :: kbound
    real, intent(in), dimension(0:4)                   :: c0
    real, intent(in), dimension(4,4)                   :: cosMoment, sinMoment
    real, intent(out), dimension(0:nx+1,0:ny+1,0:nz+1) :: phiarr
    integer                                            :: i, j, k, l, m
    real                                               :: phi, x, y, z, thet, fi, rb, costheta

    do i = 1,nx
       do j = 1,ny
          k = kbound
          phi = 0.0
          call xcoord(i,x)
          call ycoord(j,y)
          call zcoord(k,z)
          call pangle(thet,x,y,z)
          call aangle(fi,x,y)
          call distance(rb,x,y,z)
          costheta = cos(thet)
          phi = phi - c0(0)/(rb)
          do l = 1,4
             phi = phi - c0(l)*Legendre(costheta,l)/(rb**(l+1))
             do m = 1,l
                phi = phi - 2.0*(cosMoment(l,m)*cos(m*fi) + sinMoment(l,m)*sin(m*fi))*aLegendre(costheta,l,m)/(rb**(l+1))
             enddo
          enddo
          phiarr(i,j,k) = Gconst*phi
       enddo
    enddo

  end subroutine phi_zboundaries
  
#endif
  
  subroutine  poisson_solver1 (dens, phi, solver)
    
    use parameters, only: dx, dy, dz, nx, ny, nz, Gconst
    use constants, only: PI
    implicit  none
    real, parameter         ::  erro = 5e-4, gam = 1.3
    real, intent(in)        ::  dens(1:nx,1:ny,1:nz)
    real, intent(inout)     ::  phi(0:nx+1,0:ny+1,0:nz+1)
    real                    ::  Dif
    real                    ::  phi0(0:nx+1,0:ny+1,0:nz+1)
    integer                 ::  i, j, k
    integer                 ::  iter
    integer, intent(in)     :: solver
    integer, parameter      :: JACOBI = 1, GAUSSEIDEL = 2, SOR = 3
    
    iter = 0
    phi0 = phi
    Dif = 1000
    
    do while (Dif>erro)
       iter = iter +1
       Dif = 0
       
       phi0(0,:,:) = 0.0
       phi0(nx+1,:,:) = 0.0
       phi0(:,:,0) = 0.0
       phi0(:,:,nz+1) = 0.0
       phi0(:,0,:) = 0.0
       phi0(:,ny+1,:) = 0.0
       select case(solver)
       case(JACOBI)
          
          do i = 1,nx
             do j = 1,ny
                do k = 1,nz
                   phi(i,j,k) = (1.0/6.0)*(phi0(i+1,j,k) + phi0(i-1,j,k) + phi0(i,j+1,k) + phi0(i,j-1,k) + phi0(i,j,k+1) + &
                        phi0(i,j,k-1) - 4*PI*dens(i,j,k)*dx**2)
                   Dif = Dif + abs(phi(i,j,k)-phi0(i,j,k))
                enddo
             enddo
          enddo
       case(GAUSSEIDEL)
          phi(0,:,:) = phi0(0,:,:)
          phi(:,0,:) = phi0(:,0,:)
          phi(:,:,0) = phi0(:,:,0)
          phi(nx+1,:,:) = phi0(nx+1,:,:)
          phi(:,ny+1,:) = phi0(:,ny+1,:)
          phi(:,:,nz+1) = phi0(:,:,nz+1)
          
          do i = 1,nx
             do j = 1,ny
                do k = 1,nz
                   phi(i,j,k) = (1.0/6.0)*(phi0(i+1,j,k) + phi(i-1,j,k) + phi0(i,j+1,k) + phi(i,j-1,k) + phi0(i,j,k+1) + &
                        phi(i,j,k-1) - 4*PI*dens(i,j,k)*dx**2)
                   Dif = Dif + abs(phi(i,j,k)-phi0(i,j,k))
                enddo
             enddo
          enddo
       case(SOR)
          phi(0,:,:) = phi0(0,:,:)
          phi(:,0,:) = phi0(:,0,:)
          phi(:,:,0) = phi0(:,:,0)
          phi(nx+1,:,:) = phi0(nx+1,:,:)
          phi(:,ny+1,:) = phi0(:,ny+1,:)
          phi(:,:,nz+1) = phi0(:,:,nz+1)
          
          do i = 1,nx
             do j = 1,ny
                do k = 1,nz
                   phi(i,j,k) = ((1.0/6.0)*(phi0(i+1,j,k) + phi(i-1,j,k) + phi0(i,j+1,k) + phi(i,j-1,k) + phi0(i,j,k+1) + &
                        phi(i,j,k-1) - 4*PI*dens(i,j,k)*dx**2))*gam+(1-gam)*phi(i,j,k)
                   Dif = Dif + abs(phi(i,j,k)-phi0(i,j,k))
                enddo
             enddo
          enddo
       end select
       phi0 = phi
    enddo
!    write (*,*) 'Dif = ', Dif, 'iter = ', iter
  end subroutine poisson_solver1
  
  subroutine S_Solver(up,phi,S)
    
    use parameters, only: nxmin, nxmax, nymin, nymax, nzmin, nzmax, neq, nx, ny, nz, dx, dy, dz
    implicit none
    integer                 ::  i, j, k
    real, intent(in)        ::  phi(0:nx+1,0:ny+1,0:nz+1)
    real, intent(in)        ::  up(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(out)       ::  S(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real                    ::  FG(3)
    S = 0.0
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             FG = 0.0
             FG(1) = (phi(i+1,j,k) - phi(i-1,j,k))/(2.0*dx)
             FG(2) = (phi(i,j+1,k) - phi(i,j-1,k))/(2.0*dy)
             FG(3) = (phi(i,j,k+1) - phi(i,j,k-1))/(2.0*dz)
             S(1,i,j,k) = 0.0
             S(2,i,j,k) = -up(1,i,j,k)*FG(1)
             S(3,i,j,k) = -up(1,i,j,k)*FG(2)
             S(4,i,j,k) = -up(1,i,j,k)*FG(3)
             S(5,i,j,k) = -up(2,i,j,k)*up(1,i,j,k)*FG(1) - up(3,i,j,k)*up(1,i,j,k)*FG(2) - up(4,i,j,k)*up(1,i,j,k)*FG(3)
          enddo
       enddo
    enddo
  end subroutine S_Solver

end module
