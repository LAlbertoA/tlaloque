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
  
  subroutine MultiGrid(dens, phi, box_factor)

    use globals, only: Error, Residue, ret, rank, left, right, down, top, front, back, it
    use parameters, only: nx, ny, nz, lvlm, logged, logu, Gconst, xl, xr, mpix
    use constants
    
    implicit none

    real, intent(in)                     :: dens(0:nx+1,0:ny+1,0:nz+1), box_factor
    real, intent(inout)                  :: phi(0:nx+1,0:ny+1,0:nz+1)
    real, dimension(4,4)                 :: cosMoment, sinMoment
    real, dimension(0:4)                 :: c0
    integer, parameter                   :: Solve = 1, Relax = 2
    integer                              :: nxl, nyl, nzl, lv, level, iteraciones, iter
    real                                 :: Dif, dxl

    if (mod(it,10).eq.0 .and. abs(box_factor - 1.0) > 1.0e-6) then
       call Moments(c0,cosMoment,sinMoment,dens, box_factor)
       if (left<0) then
          call phi_xboundaries(0,c0,cosMoment,sinMoment,phi,box_factor)
       endif
       if (right<0) then
          call phi_xboundaries(nx+1,c0,cosMoment,sinMoment,phi,box_factor)
       endif
       if (down<0) then
          call phi_yboundaries(0,c0,cosMoment,sinMoment,phi,box_factor)
       endif
       if (top<0) then
          call phi_yboundaries(ny+1,c0,cosMoment,sinMoment,phi,box_factor)
       endif
       if (front<0) then
          call phi_zboundaries(0,c0,cosMoment,sinMoment,phi,box_factor)
       endif
       if (back<0) then
          call phi_zboundaries(nz+1,c0,cosMoment,sinMoment,phi,box_factor)
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

       dxl = box_factor*(xr-xl)/(mpix*nx)
       
       call poisson_solver(Residue(0)%data,Error(0)%data,nx,ny,nz,dxl,Relax,iteraciones,Dif,level)

       do level = 1,lvlm

          nxl = int(nx/2**level)
          nyl = int(ny/2**level)
          nzl = int(nz/2**level)

          dxl = box_factor*(xr-xl)/(mpix*nxl)

!          print*, "Restriction on ", level-1

          if (level == lvlm) then
             call resid(nxl*2,nyl*2,nzl*2,dxl/2,Error(level-1)%data,Residue(level-1)%data,ret(level-1)%data)
             call restriction(nxl*2,nyl*2,nzl*2,nxl,nyl,nzl,ret(level-1)%data,Residue(level)%data)
!             print*, "Solving residue equation on level", level
             call poisson_solver(Residue(level)%data,Error(level)%data,nxl,nyl,nzl,dxl,Solve,iteraciones,Dif,level)
          else
             call resid(nxl*2,nyl*2,nzl*2,dxl/2,Error(level-1)%data,Residue(level-1)%data,ret(level-1)%data)
             call restriction(nxl*2,nyl*2,nzl*2,nxl,nyl,nzl,ret(level-1)%data,Residue(level)%data)
!             print*, "Relaxing residue equation on level", level
             call poisson_solver(Residue(level)%data,Error(level)%data,nxl,nyl,nzl,dxl,Relax,iteraciones,Dif,level)

          endif

       enddo

       do level = lvlm-1, 0, -1

          nxl = int(nx/2**level)
          nyl = int(ny/2**level)
          nzl = int(nz/2**level)

          dxl = box_factor*(xr-xl)/(mpix*nxl)

          ret(level)%data = 0.0

          call prolongation(int(nxl/2.0),int(nyl/2.0),int(nzl/2.0),nxl,nyl,nzl,Error(level+1)%data,ret(level)%data)
!          print*, "Prolongation on ", level
          Error(level)%data = Error(level)%data + ret(level)%data
!          print*, "Relaxing residue equation on level ", level-1
          call poisson_solver(Residue(level)%data,Error(level)%data,nxl,nyl,nzl,dxl,Relax,iteraciones,Dif,level)

       enddo

       if (rank == 0 .and. logged .eqv. .true.) then     
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

  subroutine resid(nx,ny,nz,dxl,var,source,resi)

    implicit none

    integer, intent(in)   :: nx, ny, nz
    real, intent(in)      :: dxl
    real, intent(in)      :: var(0:nx+1,0:ny+1,0:nz+1), source(0:nx+1,0:ny+1,0:nz+1)
    real, intent(out)     :: resi(0:nx+1,0:ny+1,0:nz+1)
    real                  :: grad(0:nx+1,0:ny+1,0:nz+1)
    integer               :: i, j, k

    grad = 0.0

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

  subroutine  poisson_solver (source,var,nx,ny,nz,dxl,relsolv,it,Dif,level)

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
    real, intent(in)          ::  dxl
    real, parameter           ::  gam = 1.      !! gam = 1 -> Gauss-Seidel, 1 < gam < 2 -> SOR
    integer, parameter        ::  Solve = 1, Relax = 2

    real                      ::  er, d1, difp, w
    integer                   ::  i, j, k, iter, start

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

  subroutine Moments(c0,cm,sm,rho, box_factor)

    use parameters, only: nx, ny, nz, dx, dy, dz, xl, xr, yl, yr, zl, zr, mpix, mpiy, mpiz, mpi_real_kind
    use globals, only: comm3d, err, coords

    implicit none
#ifdef MPIP
    include "mpif.h"
#endif
    real, intent(in), dimension(0:nx+1,0:ny+1,0:nz+1) :: rho
    real, intent(in)                                  :: box_factor
    real, intent(out), dimension(4,4)                 :: cm, sm
    real, intent(out), dimension(0:4)                 :: c0
    real                                              :: r, coth, x, y, z, smoment
    real                                              :: c, cmp, smp, thet, fi, c0mp, cmoment, c0moment
    real                                              :: xlc, xrc, ylc, yrc, zlc, zrc
    real                                              :: dxc, dyc, dzc, dV
    real                                              :: xleft, yleft, zleft
    integer                                           :: i, j, k, l, m

    c0 = 0.0
    cm = 0.0
    sm = 0.0

    dxc = box_factor*dx
    dyc = box_factor*dy
    dzc = box_factor*dz
    dV = dxc*dyc*dzc

    xlc = 0.5*(xl+xr) - 0.5*box_factor*(xr-xl)
    xrc = 0.5*(xl+xr) + 0.5*box_factor*(xr-xl)

    ylc = 0.5*(yl+yr) - 0.5*box_factor*(yr-yl)
    yrc = 0.5*(yl+yr) + 0.5*box_factor*(yr-yl)

    zlc = 0.5*(zl+zr) - 0.5*box_factor*(zr-zl)
    zrc = 0.5*(zl+zr) + 0.5*box_factor*(zr-zl)

    xleft = xlc + coords(0)*(xrc-xlc)/mpix
    yleft = ylc + coords(1)*(yrc-ylc)/mpiy
    zleft = zlc + coords(2)*(zrc-zlc)/mpiz

    do l = 1,4
       do m = 1,l
          c0mp = 0.0
          cmp = 0.0
          smp = 0.0
          c = 2.0*Gamma(l-m+1.0)/Gamma(l+m+1.0)
          do i = 1,nx
             do j = 1,ny
                do k = 1,nz

                   x = xleft + (i - 0.5)*dxc
                   y = yleft + (j - 0.5)*dyc
                   z = zleft + (k - 0.5)*dzc
                   call pangle(thet,x,y,z)
                   call aangle(fi,x,y)
                   call distance(r,x,y,z)
                   coth = cos(thet)
                   c0mp = c0mp + Legendre(coth,l)*(r**l)*rho(i,j,k)*dV
                   cmp = cmp + aLegendre(coth,l,m)*cos(m*fi)*(r**l)*rho(i,j,k)*dV
                   smp = smp + aLegendre(coth,l,m)*sin(m*fi)*(r**l)*rho(i,j,k)*dV
                   !             print*, cmoment, smoment
                enddo
             enddo
          enddo
          cmoment = cmp
          smoment = smp
#ifdef MPIP          
          call mpi_allreduce(cmp, cmoment, 1, mpi_real_kind, mpi_sum, comm3d, err)
          call mpi_allreduce(smp, smoment, 1, mpi_real_kind, mpi_sum, comm3d, err)
#endif
          cm(l,m) = c*cmoment
          sm(l,m) = c*smoment
       enddo
#ifdef MPIP          
       call mpi_allreduce(c0mp, c0moment, 1, mpi_real_kind, mpi_sum, comm3d, err)
#else
       c0moment = c0mp
#endif
       c0(l) = c0moment
    enddo
    c0mp = 0.0
    l = 0
    m = 0
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz

             x = xleft + (i - 0.5)*dxc
             y = yleft + (j - 0.5)*dyc
             z = zleft + (k - 0.5)*dzc
             call pangle(thet,x,y,z)
             call aangle(fi,x,y)
             call distance(r,x,y,z)
             coth = cos(thet)
             c0mp = c0mp + Legendre(coth,l)*(r**l)*rho(i,j,k)*dV
          enddo
       enddo
    enddo
#ifdef MPIP
    call mpi_allreduce(c0mp, c0moment, 1, mpi_real_kind, mpi_sum, comm3d, err)
#else
    c0moment = c0mp
#endif
    c0(l) = c0moment

  end subroutine Moments

  subroutine phi_xboundaries(ibound,c0,cosMoment,sinMoment,phiarr,box_factor)

    use parameters, only: nx, ny, nz, dx, dy, dz, Gconst, xl, xr, yl, yr, zl, zr, mpix, mpiy, mpiz
    use globals, only: coords

    implicit none

    integer, intent(in)                                 :: ibound
    real, intent(in), dimension(0:4)                    :: c0
    real, intent(in), dimension(4,4)                    :: cosMoment, sinMoment
    real, intent(inout), dimension(0:nx+1,0:ny+1,0:nz+1):: phiarr
    real, intent(in)                                    :: box_factor
    real                                                :: xlc, xrc, ylc, yrc, zlc, zrc
    real                                                :: dxc, dyc, dzc
    real                                                :: xleft, yleft, zleft
    integer                                             :: i, j, k, l, m
    real                                                :: phi, x, y, z, thet, fi, rb, costheta

    dxc = box_factor*dx
    dyc = box_factor*dy
    dzc = box_factor*dz

    xlc = 0.5*(xl+xr) - 0.5*box_factor*(xr-xl)
    xrc = 0.5*(xl+xr) + 0.5*box_factor*(xr-xl)

    ylc = 0.5*(yl+yr) - 0.5*box_factor*(yr-yl)
    yrc = 0.5*(yl+yr) + 0.5*box_factor*(yr-yl)

    zlc = 0.5*(zl+zr) - 0.5*box_factor*(zr-zl)
    zrc = 0.5*(zl+zr) + 0.5*box_factor*(zr-zl)

    xleft = xlc + coords(0)*(xrc-xlc)/mpix
    yleft = ylc + coords(1)*(yrc-ylc)/mpiy
    zleft = zlc + coords(2)*(zrc-zlc)/mpiz

    i = ibound
    do j = 1,ny
       do k = 1,nz

          phi = 0.0
          x = xleft + (i - 0.5)*dxc
          y = yleft + (j - 0.5)*dyc
          z = zleft + (k - 0.5)*dzc
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

  subroutine phi_yboundaries(jbound,c0,cosMoment,sinMoment,phiarr,box_factor)

    use parameters, only: nx, ny, nz, dx, dy, dz, Gconst, xl, xr, yl, yr, zl, zr, mpix, mpiy, mpiz
    use globals, only: coords
    
    implicit none

    integer, intent(in)                                 :: jbound
    real, intent(in), dimension(0:4)                    :: c0
    real, intent(in), dimension(4,4)                    :: cosMoment, sinMoment
    real, intent(inout), dimension(0:nx+1,0:ny+1,0:nz+1):: phiarr
    real, intent(in)                                    :: box_factor
    real                                                :: xlc, xrc, ylc, yrc, zlc, zrc
    real                                                :: dxc, dyc, dzc
    real                                                :: xleft, yleft, zleft
    integer                                             :: i, j, k, l, m
    real                                                :: phi, x, y, z, thet, fi, rb, costheta

    dxc = box_factor*dx
    dyc = box_factor*dy
    dzc = box_factor*dz

    xlc = 0.5*(xl+xr) - 0.5*box_factor*(xr-xl)
    xrc = 0.5*(xl+xr) + 0.5*box_factor*(xr-xl)

    ylc = 0.5*(yl+yr) - 0.5*box_factor*(yr-yl)
    yrc = 0.5*(yl+yr) + 0.5*box_factor*(yr-yl)

    zlc = 0.5*(zl+zr) - 0.5*box_factor*(zr-zl)
    zrc = 0.5*(zl+zr) + 0.5*box_factor*(zr-zl)

    xleft = xlc + coords(0)*(xrc-xlc)/mpix
    yleft = ylc + coords(1)*(yrc-ylc)/mpiy
    zleft = zlc + coords(2)*(zrc-zlc)/mpiz

    j = jbound
    do i = 1,nx
       do k = 1,nz

          phi = 0.0
          x = xleft + (i - 0.5)*dxc
          y = yleft + (j - 0.5)*dyc
          z = zleft + (k - 0.5)*dzc
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

  subroutine phi_zboundaries(kbound,c0,cosMoment,sinMoment,phiarr,box_factor)

    use parameters, only: nx, ny, nz, dx, dy, dz, Gconst, xl, xr, yl, yr, zl, zr, mpix, mpiy, mpiz
    use globals, only: coords

    implicit none

    integer, intent(in)                                 :: kbound
    real, intent(in), dimension(0:4)                    :: c0
    real, intent(in), dimension(4,4)                    :: cosMoment, sinMoment
    real, intent(inout), dimension(0:nx+1,0:ny+1,0:nz+1):: phiarr
    real, intent(in)                                    :: box_factor
    real                                                :: xlc, xrc, ylc, yrc, zlc, zrc
    real                                                :: dxc, dyc, dzc
    real                                                :: xleft, yleft, zleft
    integer                                             :: i, j, k, l, m
    real                                                :: phi, x, y, z, thet, fi, rb, costheta

    dxc = box_factor*dx
    dyc = box_factor*dy
    dzc = box_factor*dz

    xlc = 0.5*(xl+xr) - 0.5*box_factor*(xr-xl)
    xrc = 0.5*(xl+xr) + 0.5*box_factor*(xr-xl)

    ylc = 0.5*(yl+yr) - 0.5*box_factor*(yr-yl)
    yrc = 0.5*(yl+yr) + 0.5*box_factor*(yr-yl)

    zlc = 0.5*(zl+zr) - 0.5*box_factor*(zr-zl)
    zrc = 0.5*(zl+zr) + 0.5*box_factor*(zr-zl)

    xleft = xlc + coords(0)*(xrc-xlc)/mpix
    yleft = ylc + coords(1)*(yrc-ylc)/mpiy
    zleft = zlc + coords(2)*(zrc-zlc)/mpiz

    k = kbound
    do i = 1,nx
       do j = 1,ny

          phi = 0.0
          x = xleft + (i - 0.5)*dxc
          y = yleft + (j - 0.5)*dyc
          z = zleft + (k - 0.5)*dzc
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

  subroutine  restriction_density(nx,ny,nz,nxl,nyl,nzl,var,varTemp)

    use parameters, only: neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
    implicit none

    integer, intent(in)                   :: nx, ny, nz, nxl, nyl, nzl
    real, intent(in)                      :: var(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
    real, intent(inout)                   :: varTemp(0:nx+1,0:ny+1,0:nz+1)
    integer                               :: i, j, k, i1, j1, k1, nxq, nyq, nzq

    nxq = Int(nxl/2)
    nyq = Int(nyl/2)
    nzq = Int(nzl/2)

    do i = 1,nxl
       do j = 1,nyl
          do k = 1,nzl
             i1 = i*2-1
             j1 = j*2-1
             k1 = k*2-1
             varTemp(nxq+i,nyq+j,nzq+k) = (1.0/8.0)*(var(1,i1,j1,k1)+var(1,i1+1,j1,k1)+var(1,i1,j1+1,k1)+var(1,i1,j1,k1+1)+&
                var(1,i1+1,j1+1,k1)+var(1,i1+1,j1,k1+1)+var(1,i1,j1+1,k1+1)+var(1,i1+1,j1+1,k1+1))
          end do
       end do
    end do

  end subroutine restriction_density

  subroutine prolongation_to_phi(nx,ny,nz,phi_ext,phi_hydro)

    implicit none

    integer, intent(in) :: nx, ny, nz
    real, intent(in)    :: phi_ext(0:nx+1,0:ny+1,0:nz+1)
    real, intent(inout) :: phi_hydro(0:nx+1,0:ny+1,0:nz+1)

    integer :: i, j, k
    integer :: ic, jc, kc
    integer :: ifn, jfn, kfn
    integer :: nxl, nyl, nzl
    integer :: nxq, nyq, nzq
    integer :: a, b, c
    real :: wx(0:1), wy(0:1), wz(0:1)

    nxl = nx/2 + 1
    nyl = ny/2 + 1
    nzl = nz/2 + 1

    nxq = nx/4
    nyq = ny/4
    nzq = nz/4

    do i = 1,nxl
       do j = 1,nyl
          do k = 1,nzl

             ic = nxq + i - 1
             jc = nyq + j - 1
             kc = nzq + k - 1

             ifn = 2*i - 2
             jfn = 2*j - 2
             kfn = 2*k - 2

             do a = 0,1
                if (a == 0) then
                   wx(0) = 0.75
                   wx(1) = 0.25
                else
                   wx(0) = 0.25
                   wx(1) = 0.75
                endif

                do b = 0,1
                   if (b == 0) then
                      wy(0) = 0.75
                      wy(1) = 0.25
                   else
                      wy(0) = 0.25
                      wy(1) = 0.75
                   endif

                   do c = 0,1
                      if (c == 0) then
                         wz(0) = 0.75
                         wz(1) = 0.25
                      else
                         wz(0) = 0.25
                         wz(1) = 0.75
                      endif

                      phi_hydro(ifn+a,jfn+b,kfn+c) = &
                           wx(0)*wy(0)*wz(0)*phi_ext(ic  ,jc  ,kc  ) + &
                           wx(1)*wy(0)*wz(0)*phi_ext(ic+1,jc  ,kc  ) + &
                           wx(0)*wy(1)*wz(0)*phi_ext(ic  ,jc+1,kc  ) + &
                           wx(0)*wy(0)*wz(1)*phi_ext(ic  ,jc  ,kc+1) + &
                           wx(1)*wy(1)*wz(0)*phi_ext(ic+1,jc+1,kc  ) + &
                           wx(1)*wy(0)*wz(1)*phi_ext(ic+1,jc  ,kc+1) + &
                           wx(0)*wy(1)*wz(1)*phi_ext(ic  ,jc+1,kc+1) + &
                           wx(1)*wy(1)*wz(1)*phi_ext(ic+1,jc+1,kc+1)

                   enddo
                enddo
             enddo

          enddo
       enddo
    enddo

end subroutine prolongation_to_phi

#ifdef MPIP

subroutine redistribute_density(dens)

   use parameters, only: neq, nxtot, nytot, nztot, mpix, mpiy, mpiz, nx, ny, nz, mpi_real_kind
   use globals, only: comm3d, err, coords, rank

   implicit none
   
   include "mpif.h"

   real, intent(inout) :: dens(0:nx+1,0:ny+1,0:nz+1)
   
   integer :: nxh, nyh, nzh
   integer :: ioff, joff, koff
   integer :: src_x1, src_x2, src_y1, src_y2, src_z1, src_z2
   integer :: dest_x_min, dest_x_max, dest_y_min, dest_y_max, dest_z_min, dest_z_max
   integer :: dest_x, dest_y, dest_z, dest_rank, dest_coords(3)
   integer :: dst_x1, dst_x2, dst_y1, dst_y2, dst_z1, dst_z2
   integer :: ix1, ix2, iy1, iy2, iz1, iz2
   integer :: src_i1, src_i2, src_j1, src_j2, src_k1, src_k2
   integer :: dst_i1, dst_i2, dst_j1, dst_j2, dst_k1, dst_k2

   real, allocatable :: dens_old(:,:,:), sendbuf(:,:,:), recvbuf(:,:,:)

   nxh = nx/2
   nyh = ny/2
   nzh = nz/2
   
   ioff = nxtot/4
   joff = nytot/4
   koff = nztot/4

   allocate(dens_old(0:nx+1,0:ny+1,0:nz+1))
   dens_old = dens
   dens = 0.0

   !! Notes to myself, position of the density restricted block in the global grid:
   src_x1 = ioff + coords(0)*nxh + 1
   src_x2 = ioff + (coords(0)+1)*nxh

   src_y1 = joff + coords(1)*nyh + 1
   src_y2 = joff + (coords(1)+1)*nyh

   src_z1 = koff + coords(2)*nzh + 1
   src_z2 = koff + (coords(2)+1)*nzh

   !! coordinates of the mpi rank that will receive the density block:
   dest_x_min = (src_x1 - 1)/nx
   dest_x_max = (src_x2 - 1)/nx

   dest_y_min = (src_y1 - 1)/ny
   dest_y_max = (src_y2 - 1)/ny

   dest_z_min = (src_z1 - 1)/nz
   dest_z_max = (src_z2 - 1)/nz

   do dest_x = dest_x_min, dest_x_max
      do dest_y = dest_y_min, dest_y_max
         do dest_z = dest_z_min, dest_z_max

            !! destination indices total range in the global grid:
            dst_x1 = dest_x*nx + 1
            dst_x2 = (dest_x+1)*nx

            dst_y1 = dest_y*ny + 1
            dst_y2 = (dest_y+1)*ny

            dst_z1 = dest_z*nz + 1
            dst_z2 = (dest_z+1)*nz

            !! overlap of the source block and the destination block in the global grid:
            ix1 = max(src_x1,dst_x1)
            ix2 = min(src_x2,dst_x2)

            iy1 = max(src_y1,dst_y1)
            iy2 = min(src_y2,dst_y2)

            iz1 = max(src_z1,dst_z1)
            iz2 = min(src_z2,dst_z2)

            if (ix1 <= ix2 .and. iy1 <= iy2 .and. iz1 <= iz2) then
               
               !! indices of the source block in the local grid of this rank:
               src_i1 = nx/4 + 1 + (ix1 - src_x1)
               src_i2 = nx/4 + 1 + (ix2 - src_x1)

               src_j1 = ny/4 + 1 + (iy1 - src_y1)
               src_j2 = ny/4 + 1 + (iy2 - src_y1)

               src_k1 = nz/4 + 1 + (iz1 - src_z1)
               src_k2 = nz/4 + 1 + (iz2 - src_z1)

               !! indices of the destination block in the local grid of the destination rank:
               dst_i1 = ix1 - dest_x*nx
               dst_i2 = ix2 - dest_x*nx

               dst_j1 = iy1 - dest_y*ny
               dst_j2 = iy2 - dest_y*ny

               dst_k1 = iz1 - dest_z*nz
               dst_k2 = iz2 - dest_z*nz

               dest_coords = (/dest_x, dest_y, dest_z/)

               call mpi_cart_rank(comm3d, dest_coords, dest_rank, err)

               if (dest_rank == rank) then
                  dens(dst_i1:dst_i2,dst_j1:dst_j2,dst_k1:dst_k2) = &
                     dens_old(src_i1:src_i2,src_j1:src_j2,src_k1:src_k2)
               else
                  call mpi_send(dens_old(src_i1:src_i2,src_j1:src_j2,src_k1:src_k2), &
                        (src_i2-src_i1+1)*(src_j2-src_j1+1)*(src_k2-src_k1+1), mpi_real_kind, dest_rank, 0, comm3d, err)
               endif

               endif
            enddo
         enddo
      enddo
      
      deallocate (dens_old)

   end subroutine redistribute_density

subroutine setup_density_redistribution()

    use parameters, only: nx, ny, nz, mpix, mpiy, mpiz
    use globals, only: coords, rank, comm3d, err
    use density_redistribution_data

    implicit none

#ifdef MPIP

	include "mpif.h"

    integer :: NXG, NYG, NZG
    integer :: nxh, nyh, nzh
    integer :: ioff, joff, koff

    integer :: src_x1, src_x2, src_y1, src_y2, src_z1, src_z2
    integer :: dest_x_min, dest_x_max
    integer :: dest_y_min, dest_y_max
    integer :: dest_z_min, dest_z_max

    integer :: dest_x, dest_y, dest_z
    integer :: dest_rank
    integer :: dest_coords(3)

    integer :: dst_x1, dst_x2, dst_y1, dst_y2, dst_z1, dst_z2
    integer :: ix1, ix2, iy1, iy2, iz1, iz2

    integer :: src_i1, src_i2, src_j1, src_j2, src_k1, src_k2
    integer :: dst_i1, dst_i2, dst_j1, dst_j2, dst_k1, dst_k2

    integer :: nprocs
    integer :: b, p, offset
    integer :: nblocks_p
    integer :: recv_count

    integer :: my_meta(1 + meta_size*max_send_blocks)
    integer, allocatable :: all_meta(:,:)

    integer :: tmp_send_dest(max_send_blocks)
    integer :: tmp_send_nblock(max_send_blocks)

    integer :: tmp_send_src_i1(max_send_blocks), tmp_send_src_i2(max_send_blocks)
    integer :: tmp_send_src_j1(max_send_blocks), tmp_send_src_j2(max_send_blocks)
    integer :: tmp_send_src_k1(max_send_blocks), tmp_send_src_k2(max_send_blocks)

    integer :: tmp_send_dst_i1(max_send_blocks), tmp_send_dst_i2(max_send_blocks)
    integer :: tmp_send_dst_j1(max_send_blocks), tmp_send_dst_j2(max_send_blocks)
    integer :: tmp_send_dst_k1(max_send_blocks), tmp_send_dst_k2(max_send_blocks)

    if (mod(nx,4) /= 0 .or. mod(ny,4) /= 0 .or. mod(nz,4) /= 0) then
       if (rank == 0) then
          print *, "ERROR: setup_density_redistribution requires nx,ny,nz divisible by 4"
       endif
       stop
    endif

    call MPI_Comm_size(comm3d, nprocs, err)

    NXG = nx*mpix
    NYG = ny*mpiy
    NZG = nz*mpiz

    if (mod(NXG,4) /= 0 .or. mod(NYG,4) /= 0 .or. mod(NZG,4) /= 0) then
       if (rank == 0) then
          print *, "ERROR: global grid must be divisible by 4 in each direction"
       endif
       stop
    endif

    nxh = nx/2
    nyh = ny/2
    nzh = nz/2

    ioff = NXG/4
    joff = NYG/4
    koff = NZG/4

    ! ============================================================
    ! 1. Calculate local send metadata
    ! ============================================================

    nsend_rho_blocks = 0

    src_x1 = ioff + coords(0)*nxh + 1
    src_x2 = ioff + (coords(0)+1)*nxh

    src_y1 = joff + coords(1)*nyh + 1
    src_y2 = joff + (coords(1)+1)*nyh

    src_z1 = koff + coords(2)*nzh + 1
    src_z2 = koff + (coords(2)+1)*nzh

    dest_x_min = (src_x1 - 1)/nx
    dest_x_max = (src_x2 - 1)/nx

    dest_y_min = (src_y1 - 1)/ny
    dest_y_max = (src_y2 - 1)/ny

    dest_z_min = (src_z1 - 1)/nz
    dest_z_max = (src_z2 - 1)/nz

    do dest_x = dest_x_min, dest_x_max
       do dest_y = dest_y_min, dest_y_max
          do dest_z = dest_z_min, dest_z_max

             dst_x1 = dest_x*nx + 1
             dst_x2 = (dest_x+1)*nx

             dst_y1 = dest_y*ny + 1
             dst_y2 = (dest_y+1)*ny

             dst_z1 = dest_z*nz + 1
             dst_z2 = (dest_z+1)*nz

             ix1 = max(src_x1,dst_x1)
             ix2 = min(src_x2,dst_x2)

             iy1 = max(src_y1,dst_y1)
             iy2 = min(src_y2,dst_y2)

             iz1 = max(src_z1,dst_z1)
             iz2 = min(src_z2,dst_z2)

             if (ix1 <= ix2 .and. iy1 <= iy2 .and. iz1 <= iz2) then

                nsend_rho_blocks = nsend_rho_blocks + 1

                if (nsend_rho_blocks > max_send_blocks) then
                   print *, "ERROR: too many send blocks in setup_density_redistribution"
                   stop
                endif

                src_i1 = nx/4 + 1 + (ix1 - src_x1)
                src_i2 = nx/4 + 1 + (ix2 - src_x1)

                src_j1 = ny/4 + 1 + (iy1 - src_y1)
                src_j2 = ny/4 + 1 + (iy2 - src_y1)

                src_k1 = nz/4 + 1 + (iz1 - src_z1)
                src_k2 = nz/4 + 1 + (iz2 - src_z1)

                dst_i1 = ix1 - dest_x*nx
                dst_i2 = ix2 - dest_x*nx

                dst_j1 = iy1 - dest_y*ny
                dst_j2 = iy2 - dest_y*ny

                dst_k1 = iz1 - dest_z*nz
                dst_k2 = iz2 - dest_z*nz

                dest_coords = (/dest_x,dest_y,dest_z/)
                call MPI_Cart_rank(comm3d,dest_coords,dest_rank,err)

                tmp_send_dest(nsend_rho_blocks) = dest_rank

                tmp_send_src_i1(nsend_rho_blocks) = src_i1
                tmp_send_src_i2(nsend_rho_blocks) = src_i2
                tmp_send_src_j1(nsend_rho_blocks) = src_j1
                tmp_send_src_j2(nsend_rho_blocks) = src_j2
                tmp_send_src_k1(nsend_rho_blocks) = src_k1
                tmp_send_src_k2(nsend_rho_blocks) = src_k2

                tmp_send_dst_i1(nsend_rho_blocks) = dst_i1
                tmp_send_dst_i2(nsend_rho_blocks) = dst_i2
                tmp_send_dst_j1(nsend_rho_blocks) = dst_j1
                tmp_send_dst_j2(nsend_rho_blocks) = dst_j2
                tmp_send_dst_k1(nsend_rho_blocks) = dst_k1
                tmp_send_dst_k2(nsend_rho_blocks) = dst_k2

                tmp_send_nblock(nsend_rho_blocks) = &
                     (src_i2-src_i1+1) * &
                     (src_j2-src_j1+1) * &
                     (src_k2-src_k1+1)

             endif

          enddo
       enddo
    enddo

    ! Allocate exact-size send metadata arrays.

    allocate(send_dest(nsend_rho_blocks))
    allocate(send_block_id(nsend_rho_blocks))
    allocate(send_nblock(nsend_rho_blocks))

    allocate(send_src_i1(nsend_rho_blocks), send_src_i2(nsend_rho_blocks))
    allocate(send_src_j1(nsend_rho_blocks), send_src_j2(nsend_rho_blocks))
    allocate(send_src_k1(nsend_rho_blocks), send_src_k2(nsend_rho_blocks))

    allocate(send_dst_i1(nsend_rho_blocks), send_dst_i2(nsend_rho_blocks))
    allocate(send_dst_j1(nsend_rho_blocks), send_dst_j2(nsend_rho_blocks))
    allocate(send_dst_k1(nsend_rho_blocks), send_dst_k2(nsend_rho_blocks))

    do b = 1, nsend_rho_blocks
       send_dest(b) = tmp_send_dest(b)
       send_block_id(b) = b
       send_nblock(b) = tmp_send_nblock(b)

       send_src_i1(b) = tmp_send_src_i1(b)
       send_src_i2(b) = tmp_send_src_i2(b)
       send_src_j1(b) = tmp_send_src_j1(b)
       send_src_j2(b) = tmp_send_src_j2(b)
       send_src_k1(b) = tmp_send_src_k1(b)
       send_src_k2(b) = tmp_send_src_k2(b)

       send_dst_i1(b) = tmp_send_dst_i1(b)
       send_dst_i2(b) = tmp_send_dst_i2(b)
       send_dst_j1(b) = tmp_send_dst_j1(b)
       send_dst_j2(b) = tmp_send_dst_j2(b)
       send_dst_k1(b) = tmp_send_dst_k1(b)
       send_dst_k2(b) = tmp_send_dst_k2(b)
    enddo

    ! ============================================================
    ! 2. Share metadata with all ranks
    ! ============================================================

    my_meta = 0
    my_meta(1) = nsend_rho_blocks

    do b = 1, nsend_rho_blocks

       offset = 1 + (b-1)*meta_size

       my_meta(offset+1)  = send_dest(b)
       my_meta(offset+2)  = send_src_i1(b)
       my_meta(offset+3)  = send_src_i2(b)
       my_meta(offset+4)  = send_src_j1(b)
       my_meta(offset+5)  = send_src_j2(b)
       my_meta(offset+6)  = send_src_k1(b)
       my_meta(offset+7)  = send_src_k2(b)
       my_meta(offset+8)  = send_dst_i1(b)
       my_meta(offset+9)  = send_dst_i2(b)
       my_meta(offset+10) = send_dst_j1(b)
       my_meta(offset+11) = send_dst_j2(b)
       my_meta(offset+12) = send_dst_k1(b)
       my_meta(offset+13) = send_dst_k2(b)
       my_meta(offset+14) = send_nblock(b)

    enddo

    allocate(all_meta(1 + meta_size*max_send_blocks, nprocs))

    call MPI_Allgather(my_meta, 1 + meta_size*max_send_blocks, MPI_INTEGER, &
                       all_meta, 1 + meta_size*max_send_blocks, MPI_INTEGER, &
                       comm3d, err)

    ! ============================================================
    ! 3. Count receives for this rank
    ! ============================================================

    recv_count = 0

    do p = 1, nprocs
       nblocks_p = all_meta(1,p)

       do b = 1, nblocks_p
          offset = 1 + (b-1)*meta_size

          if (all_meta(offset+1,p) == rank .and. (p-1) /= rank) then
             recv_count = recv_count + 1
          endif
       enddo
    enddo

    nrecv_rho_blocks = recv_count

    allocate(recv_source(nrecv_rho_blocks))
    allocate(recv_block_id(nrecv_rho_blocks))
    allocate(recv_nblock(nrecv_rho_blocks))

    allocate(recv_dst_i1(nrecv_rho_blocks), recv_dst_i2(nrecv_rho_blocks))
    allocate(recv_dst_j1(nrecv_rho_blocks), recv_dst_j2(nrecv_rho_blocks))
    allocate(recv_dst_k1(nrecv_rho_blocks), recv_dst_k2(nrecv_rho_blocks))

    recv_count = 0

    do p = 1, nprocs
       nblocks_p = all_meta(1,p)

       do b = 1, nblocks_p
          offset = 1 + (b-1)*meta_size

          if (all_meta(offset+1,p) == rank .and. (p-1) /= rank) then

             recv_count = recv_count + 1

             recv_source(recv_count)   = p - 1
             recv_block_id(recv_count) = b
             recv_nblock(recv_count)   = all_meta(offset+14,p)

             recv_dst_i1(recv_count) = all_meta(offset+8,p)
             recv_dst_i2(recv_count) = all_meta(offset+9,p)
             recv_dst_j1(recv_count) = all_meta(offset+10,p)
             recv_dst_j2(recv_count) = all_meta(offset+11,p)
             recv_dst_k1(recv_count) = all_meta(offset+12,p)
             recv_dst_k2(recv_count) = all_meta(offset+13,p)

          endif
       enddo
    enddo

    ! ============================================================
    ! 4. Allocate persistent buffers and requests
    ! ============================================================

    allocate(sendbuf(nsend_rho_blocks))
    allocate(recvbuf(nrecv_rho_blocks))

    do b = 1, nsend_rho_blocks
       if (send_dest(b) /= rank) then
          allocate(sendbuf(b)%data(send_nblock(b)))
       endif
    enddo

    do b = 1, nrecv_rho_blocks
       allocate(recvbuf(b)%data(recv_nblock(b)))
    enddo

    allocate(send_req(nsend_rho_blocks))
    allocate(recv_req(nrecv_rho_blocks))

    density_redistribution_ready = .true.

    deallocate(all_meta)

#endif

end subroutine setup_density_redistribution

subroutine setup_phi_prolongation_redistribution()

    use parameters, only: nx, ny, nz, mpix, mpiy, mpiz
    use globals, only: coords, rank, comm3d, err
    use density_redistribution_data

    implicit none

#ifdef MPIP

	include "mpif.h"

    integer :: NXG, NYG, NZG
    integer :: nxh, nyh, nzh
    integer :: nxq, nyq, nzq
    integer :: ioff, joff, koff

    integer :: need_x1, need_x2
    integer :: need_y1, need_y2
    integer :: need_z1, need_z2

    integer :: src_x_min, src_x_max
    integer :: src_y_min, src_y_max
    integer :: src_z_min, src_z_max

    integer :: src_x, src_y, src_z
    integer :: src_rank
    integer :: src_coords(3)

    integer :: owner_x1, owner_x2
    integer :: owner_y1, owner_y2
    integer :: owner_z1, owner_z2

    integer :: ix1, ix2, iy1, iy2, iz1, iz2

    integer :: src_i1, src_i2, src_j1, src_j2, src_k1, src_k2
    integer :: dst_i1, dst_i2, dst_j1, dst_j2, dst_k1, dst_k2

    integer :: nprocs
    integer :: b, p, offset, nblocks_p
    integer :: send_count

    integer :: my_meta(1 + phi_meta_size*max_phi_recv_blocks)
    integer, allocatable :: all_meta(:,:)

    call MPI_Comm_size(comm3d, nprocs, err)

    NXG = nx*mpix
    NYG = ny*mpiy
    NZG = nz*mpiz

    nxh = nx/2
    nyh = ny/2
    nzh = nz/2

    nxq = nx/4
    nyq = ny/4
    nzq = nz/4

    ioff = NXG/4
    joff = NYG/4
    koff = NZG/4

    ! ============================================================
    ! 1. Este rank hidro calcula qué bloques de PHI_ext necesita.
    !
    ! Para la prolongación corregida que usa ic e ic+1,
    ! necesitamos llenar localmente:
    !
    !   nx/4 : nx/4 + nx/2 + 1
    !
    ! Por eso el bloque global requerido es:
    !
    !   ioff + coords*nxh : ioff + (coords+1)*nxh + 1
    ! ============================================================

    need_x1 = ioff + coords(0)*nxh
    need_x2 = ioff + (coords(0)+1)*nxh + 1

    need_y1 = joff + coords(1)*nyh
    need_y2 = joff + (coords(1)+1)*nyh + 1

    need_z1 = koff + coords(2)*nzh
    need_z2 = koff + (coords(2)+1)*nzh + 1

    src_x_min = (need_x1 - 1)/nx
    src_x_max = (need_x2 - 1)/nx

    src_y_min = (need_y1 - 1)/ny
    src_y_max = (need_y2 - 1)/ny

    src_z_min = (need_z1 - 1)/nz
    src_z_max = (need_z2 - 1)/nz

    nrecv_phi_blocks = 0

    do src_x = src_x_min, src_x_max
       do src_y = src_y_min, src_y_max
          do src_z = src_z_min, src_z_max

             owner_x1 = src_x*nx + 1
             owner_x2 = (src_x+1)*nx

             owner_y1 = src_y*ny + 1
             owner_y2 = (src_y+1)*ny

             owner_z1 = src_z*nz + 1
             owner_z2 = (src_z+1)*nz

             ix1 = max(need_x1, owner_x1)
             ix2 = min(need_x2, owner_x2)

             iy1 = max(need_y1, owner_y1)
             iy2 = min(need_y2, owner_y2)

             iz1 = max(need_z1, owner_z1)
             iz2 = min(need_z2, owner_z2)

             if (ix1 <= ix2 .and. iy1 <= iy2 .and. iz1 <= iz2) then

                nrecv_phi_blocks = nrecv_phi_blocks + 1

                if (nrecv_phi_blocks > max_phi_recv_blocks) then
                   print *, "ERROR: too many phi recv blocks"
                   stop
                endif

             endif

          enddo
       enddo
    enddo

    allocate(phi_recv_source(nrecv_phi_blocks))
    allocate(phi_recv_block_id(nrecv_phi_blocks))
    allocate(phi_recv_nblock(nrecv_phi_blocks))

    allocate(phi_recv_dst_i1(nrecv_phi_blocks), phi_recv_dst_i2(nrecv_phi_blocks))
    allocate(phi_recv_dst_j1(nrecv_phi_blocks), phi_recv_dst_j2(nrecv_phi_blocks))
    allocate(phi_recv_dst_k1(nrecv_phi_blocks), phi_recv_dst_k2(nrecv_phi_blocks))

	allocate(phi_recv_src_i1(nrecv_phi_blocks), phi_recv_src_i2(nrecv_phi_blocks))
	allocate(phi_recv_src_j1(nrecv_phi_blocks), phi_recv_src_j2(nrecv_phi_blocks))
	allocate(phi_recv_src_k1(nrecv_phi_blocks), phi_recv_src_k2(nrecv_phi_blocks))

    my_meta = 0
    my_meta(1) = nrecv_phi_blocks

    b = 0

    do src_x = src_x_min, src_x_max
       do src_y = src_y_min, src_y_max
          do src_z = src_z_min, src_z_max

             owner_x1 = src_x*nx + 1
             owner_x2 = (src_x+1)*nx

             owner_y1 = src_y*ny + 1
             owner_y2 = (src_y+1)*ny

             owner_z1 = src_z*nz + 1
             owner_z2 = (src_z+1)*nz

             ix1 = max(need_x1, owner_x1)
             ix2 = min(need_x2, owner_x2)

             iy1 = max(need_y1, owner_y1)
             iy2 = min(need_y2, owner_y2)

             iz1 = max(need_z1, owner_z1)
             iz2 = min(need_z2, owner_z2)

             if (ix1 <= ix2 .and. iy1 <= iy2 .and. iz1 <= iz2) then

                b = b + 1

                src_coords = (/src_x,src_y,src_z/)
                call MPI_Cart_rank(comm3d, src_coords, src_rank, err)

                ! Índices fuente en el rank que posee PHI extendido.
                src_i1 = ix1 - src_x*nx
                src_i2 = ix2 - src_x*nx

                src_j1 = iy1 - src_y*ny
                src_j2 = iy2 - src_y*ny

                src_k1 = iz1 - src_z*nz
                src_k2 = iz2 - src_z*nz

                ! Índices destino en phi_ext_local del rank hidro.
                dst_i1 = nxq + (ix1 - need_x1)
                dst_i2 = nxq + (ix2 - need_x1)

                dst_j1 = nyq + (iy1 - need_y1)
                dst_j2 = nyq + (iy2 - need_y1)

                dst_k1 = nzq + (iz1 - need_z1)
                dst_k2 = nzq + (iz2 - need_z1)

                phi_recv_source(b) = src_rank
                phi_recv_block_id(b) = b

                phi_recv_dst_i1(b) = dst_i1
                phi_recv_dst_i2(b) = dst_i2
                phi_recv_dst_j1(b) = dst_j1
                phi_recv_dst_j2(b) = dst_j2
                phi_recv_dst_k1(b) = dst_k1
                phi_recv_dst_k2(b) = dst_k2

				phi_recv_src_i1(b) = src_i1
				phi_recv_src_i2(b) = src_i2

				phi_recv_src_j1(b) = src_j1
				phi_recv_src_j2(b) = src_j2

				phi_recv_src_k1(b) = src_k1
				phi_recv_src_k2(b) = src_k2

                phi_recv_nblock(b) = &
                     (src_i2-src_i1+1) * &
                     (src_j2-src_j1+1) * &
                     (src_k2-src_k1+1)

                ! Metadata que todos compartirán.
                ! Ojo: desde el punto de vista global,
                ! esta metadata describe una petición:
                !
                !   source rank = src_rank
                !   destination rank = rank
                !
                offset = 1 + (b-1)*phi_meta_size

                my_meta(offset+1)  = src_rank
                my_meta(offset+2)  = rank

                my_meta(offset+3)  = src_i1
                my_meta(offset+4)  = src_i2
                my_meta(offset+5)  = src_j1
                my_meta(offset+6)  = src_j2
                my_meta(offset+7)  = src_k1
                my_meta(offset+8)  = src_k2

                my_meta(offset+9)  = dst_i1
                my_meta(offset+10) = dst_i2
                my_meta(offset+11) = dst_j1
                my_meta(offset+12) = dst_j2
                my_meta(offset+13) = dst_k1
                my_meta(offset+14) = dst_k2

             endif

          enddo
       enddo
    enddo

    ! ============================================================
    ! 2. Compartir las peticiones con todos los ranks.
    ! ============================================================

    allocate(all_meta(1 + phi_meta_size*max_phi_recv_blocks, nprocs))

    call MPI_Allgather(my_meta, 1 + phi_meta_size*max_phi_recv_blocks, MPI_INTEGER, &
                       all_meta, 1 + phi_meta_size*max_phi_recv_blocks, MPI_INTEGER, &
                       comm3d, err)

    ! ============================================================
    ! 3. Contar cuántos bloques debe enviar este rank.
    ! ============================================================

    send_count = 0

    do p = 1, nprocs

       nblocks_p = all_meta(1,p)

       do b = 1, nblocks_p

          offset = 1 + (b-1)*phi_meta_size

          src_rank = all_meta(offset+1,p)

          if (src_rank == rank .and. (p-1) /= rank) then
             send_count = send_count + 1
          endif

       enddo
    enddo

    nsend_phi_blocks = send_count

    allocate(phi_send_dest(nsend_phi_blocks))
    allocate(phi_send_block_id(nsend_phi_blocks))
    allocate(phi_send_nblock(nsend_phi_blocks))

    allocate(phi_send_src_i1(nsend_phi_blocks), phi_send_src_i2(nsend_phi_blocks))
    allocate(phi_send_src_j1(nsend_phi_blocks), phi_send_src_j2(nsend_phi_blocks))
    allocate(phi_send_src_k1(nsend_phi_blocks), phi_send_src_k2(nsend_phi_blocks))

    allocate(phi_send_dst_i1(nsend_phi_blocks), phi_send_dst_i2(nsend_phi_blocks))
    allocate(phi_send_dst_j1(nsend_phi_blocks), phi_send_dst_j2(nsend_phi_blocks))
    allocate(phi_send_dst_k1(nsend_phi_blocks), phi_send_dst_k2(nsend_phi_blocks))

    send_count = 0

    do p = 1, nprocs

       nblocks_p = all_meta(1,p)

       do b = 1, nblocks_p

          offset = 1 + (b-1)*phi_meta_size

          src_rank = all_meta(offset+1,p)

          if (src_rank == rank .and. (p-1) /= rank) then

             send_count = send_count + 1

             phi_send_dest(send_count) = all_meta(offset+2,p)
             phi_send_block_id(send_count) = b

             phi_send_src_i1(send_count) = all_meta(offset+3,p)
             phi_send_src_i2(send_count) = all_meta(offset+4,p)
             phi_send_src_j1(send_count) = all_meta(offset+5,p)
             phi_send_src_j2(send_count) = all_meta(offset+6,p)
             phi_send_src_k1(send_count) = all_meta(offset+7,p)
             phi_send_src_k2(send_count) = all_meta(offset+8,p)

             phi_send_dst_i1(send_count) = all_meta(offset+9,p)
             phi_send_dst_i2(send_count) = all_meta(offset+10,p)
             phi_send_dst_j1(send_count) = all_meta(offset+11,p)
             phi_send_dst_j2(send_count) = all_meta(offset+12,p)
             phi_send_dst_k1(send_count) = all_meta(offset+13,p)
             phi_send_dst_k2(send_count) = all_meta(offset+14,p)

             phi_send_nblock(send_count) = &
                  (phi_send_src_i2(send_count)-phi_send_src_i1(send_count)+1) * &
                  (phi_send_src_j2(send_count)-phi_send_src_j1(send_count)+1) * &
                  (phi_send_src_k2(send_count)-phi_send_src_k1(send_count)+1)

          endif

       enddo
    enddo

    ! ============================================================
    ! 4. Buffers persistentes.
    ! ============================================================

    allocate(phi_recvbuf(nrecv_phi_blocks))
    allocate(phi_sendbuf(nsend_phi_blocks))

    do b = 1, nrecv_phi_blocks
       if (phi_recv_source(b) /= rank) then
          allocate(phi_recvbuf(b)%data(phi_recv_nblock(b)))
       endif
    enddo

    do b = 1, nsend_phi_blocks
       allocate(phi_sendbuf(b)%data(phi_send_nblock(b)))
    enddo

    allocate(phi_recv_req(nrecv_phi_blocks))
    allocate(phi_send_req(nsend_phi_blocks))

    phi_redistribution_ready = .true.

    deallocate(all_meta)

#endif

end subroutine setup_phi_prolongation_redistribution

subroutine redistribute_density_blocks(dens)

    use parameters, only: nx, ny, nz, mpi_real_kind
    use globals, only: rank, comm3d, err
    use density_redistribution_data

    implicit none

    real, intent(inout) :: dens(0:nx+1,0:ny+1,0:nz+1)

#ifdef MPIP

	include "mpif.h"

    real, allocatable :: dens_old(:,:,:)
    integer :: b, i, j, k, n
    integer :: tag

    if (.not. density_redistribution_ready) then
       if (rank == 0) then
          print *, "ERROR: redistribute_density_blocks called before setup_density_redistribution"
       endif
       stop
    endif

    allocate(dens_old(0:nx+1,0:ny+1,0:nz+1))

    dens_old = dens
    dens = 0.0

    ! ============================================================
    ! 1. Post receives
    ! ============================================================

    do b = 1, nrecv_rho_blocks

       tag = tag_base + recv_block_id(b)

       call MPI_Irecv(recvbuf(b)%data, recv_nblock(b), mpi_real_kind, &
                      recv_source(b), tag, comm3d, recv_req(b), err)

    enddo

    ! ============================================================
    ! 2. Local copies and nonlocal sends
    ! ============================================================

    do b = 1, nsend_rho_blocks

       if (send_dest(b) == rank) then

          dens(send_dst_i1(b):send_dst_i2(b), &
               send_dst_j1(b):send_dst_j2(b), &
               send_dst_k1(b):send_dst_k2(b)) = &
          dens_old(send_src_i1(b):send_src_i2(b), &
                   send_src_j1(b):send_src_j2(b), &
                   send_src_k1(b):send_src_k2(b))

          send_req(b) = MPI_REQUEST_NULL

       else

          n = 0

          do k = send_src_k1(b), send_src_k2(b)
             do j = send_src_j1(b), send_src_j2(b)
                do i = send_src_i1(b), send_src_i2(b)
                   n = n + 1
                   sendbuf(b)%data(n) = dens_old(i,j,k)
                enddo
             enddo
          enddo

          tag = tag_base + send_block_id(b)

          call MPI_Isend(sendbuf(b)%data, send_nblock(b), mpi_real_kind, &
                         send_dest(b), tag, comm3d, send_req(b), err)

       endif

    enddo

    ! ============================================================
    ! 3. Wait receives and unpack
    ! ============================================================

    if (nrecv_rho_blocks > 0) then
       call MPI_Waitall(nrecv_rho_blocks, recv_req, MPI_STATUSES_IGNORE, err)
    endif

    do b = 1, nrecv_rho_blocks

       n = 0

       do k = recv_dst_k1(b), recv_dst_k2(b)
          do j = recv_dst_j1(b), recv_dst_j2(b)
             do i = recv_dst_i1(b), recv_dst_i2(b)
                n = n + 1
                dens(i,j,k) = recvbuf(b)%data(n)
             enddo
          enddo
       enddo

    enddo

    ! ============================================================
    ! 4. Wait sends
    ! ============================================================

    if (nsend_rho_blocks > 0) then
       call MPI_Waitall(nsend_rho_blocks, send_req, MPI_STATUSES_IGNORE, err)
    endif

    deallocate(dens_old)

#else

    return

#endif

end subroutine redistribute_density_blocks

subroutine redistribute_phi_for_prolongation(phi_ext, phi_local)

    use parameters, only: nx, ny, nz, mpi_real_kind
    use globals, only: rank, comm3d, err
    use density_redistribution_data

    implicit none

#ifdef MPIP
    include "mpif.h"
#endif

    real, intent(in)    :: phi_ext(0:nx+1,0:ny+1,0:nz+1)
    real, intent(inout) :: phi_local(0:nx+1,0:ny+1,0:nz+1)

#ifdef MPIP

    integer :: b, i, j, k, n
    integer :: tag

    if (.not. phi_redistribution_ready) then
       if (rank == 0) then
          print *, "ERROR: redistribute_phi_for_prolongation called before setup"
       endif
       stop
    endif

    phi_local = 0.0

    ! ============================================================
    ! 1. Post receives.
    ! ============================================================

    do b = 1, nrecv_phi_blocks

       if (phi_recv_source(b) /= rank) then

          tag = phi_tag_base + phi_recv_block_id(b)

          call MPI_Irecv(phi_recvbuf(b)%data, phi_recv_nblock(b), mpi_real_kind, &
                         phi_recv_source(b), tag, comm3d, phi_recv_req(b), err)

       else

          phi_recv_req(b) = MPI_REQUEST_NULL

       endif

    enddo

    ! ============================================================
    ! 2. Local copies and sends.
    ! ============================================================

    do b = 1, nrecv_phi_blocks

          ! Local copy: este rank ya tiene el pedazo de PHI_ext que necesita.
		if (phi_recv_source(b) == rank) then

		phi_local(phi_recv_dst_i1(b):phi_recv_dst_i2(b), &
					phi_recv_dst_j1(b):phi_recv_dst_j2(b), &
					phi_recv_dst_k1(b):phi_recv_dst_k2(b)) = &
		phi_ext(phi_recv_src_i1(b):phi_recv_src_i2(b), &
				phi_recv_src_j1(b):phi_recv_src_j2(b), &
				phi_recv_src_k1(b):phi_recv_src_k2(b))

		endif

    enddo

    do b = 1, nsend_phi_blocks

       n = 0

       do k = phi_send_src_k1(b), phi_send_src_k2(b)
          do j = phi_send_src_j1(b), phi_send_src_j2(b)
             do i = phi_send_src_i1(b), phi_send_src_i2(b)
                n = n + 1
                phi_sendbuf(b)%data(n) = phi_ext(i,j,k)
             enddo
          enddo
       enddo

       tag = phi_tag_base + phi_send_block_id(b)

       call MPI_Isend(phi_sendbuf(b)%data, phi_send_nblock(b), mpi_real_kind, &
                      phi_send_dest(b), tag, comm3d, phi_send_req(b), err)

    enddo

    ! ============================================================
    ! 3. Wait receives and unpack.
    ! ============================================================

    if (nrecv_phi_blocks > 0) then
       call MPI_Waitall(nrecv_phi_blocks, phi_recv_req, MPI_STATUSES_IGNORE, err)
    endif

    do b = 1, nrecv_phi_blocks

       if (phi_recv_source(b) /= rank) then

          n = 0

          do k = phi_recv_dst_k1(b), phi_recv_dst_k2(b)
             do j = phi_recv_dst_j1(b), phi_recv_dst_j2(b)
                do i = phi_recv_dst_i1(b), phi_recv_dst_i2(b)
                   n = n + 1
                   phi_local(i,j,k) = phi_recvbuf(b)%data(n)
                enddo
             enddo
          enddo

       endif

    enddo

    ! ============================================================
    ! 4. Wait sends.
    ! ============================================================

    if (nsend_phi_blocks > 0) then
       call MPI_Waitall(nsend_phi_blocks, phi_send_req, MPI_STATUSES_IGNORE, err)
    endif

#else

    phi_local = phi_ext

#endif

end subroutine redistribute_phi_for_prolongation

#endif
end module
