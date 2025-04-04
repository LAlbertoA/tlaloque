  subroutine borders(U,order)

    use parameters, only: neq, nxmin, nymin, nzmin, nxmax, nymax, nzmax
    use user
    
    implicit none

    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(inout) :: U
    integer, intent(in)                                                        :: order

    select case(order)
    case(1)

       call boundariesI(U)

       call user_boundaries(U)
       
    case(2)

       call boundariesII(U)

       call user_boundaries(U)
       
    end select

  end subroutine borders
  
  subroutine boundariesI(U)
    
    use parameters, only: nx, ny, nz, nxmin, nymin, nzmin, nxmax, nymax, nzmax, nghost, neq, np, &
         LB, RB, TB, DB, BB, FB
#ifdef MPIP
    use parameters, only: mpi_real_kind
    use globals, only: rank, err, left, right, down, top, front, back, comm3d
#endif
    use constants
    
    implicit none

#ifdef MPIP

    include "mpif.h"

#endif
    
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(inout) :: U
    integer, parameter                      ::  nyg1 = ny+1, nzg1 = nz+1, nxg1 = nx+1

#ifdef MPIP
    
    integer                                 ::  status(MPI_STATUS_SIZE)
    real, dimension(neq,1,0:ny+1,0:nz+1)    ::  sendr,recvr,sendl,recvl
    real, dimension(neq,0:nx+1,1,0:nz+1)    ::  sendt,recvt,sendd,recvd
    real, dimension(neq,0:nx+1,0:ny+1,1)    ::  sendb,recvb,sendf,recvf
    integer, parameter                      ::  srsizex = neq*(ny+nghost)*(nz+nghost)
    integer, parameter                      ::  srsizey = neq*(nx+nghost)*(nz+nghost)
    integer, parameter                      ::  srsizez = neq*(nx+nghost)*(ny+nghost)
    
    sendr(:,1,:,:)=U(:,    nx,0:nyg1,0:nzg1)
    sendl(:,1,:,:)=U(:,     1,0:nyg1,0:nzg1)
    sendt(:,:,1,:)=U(:,0:nxg1,    ny,0:nzg1)
    sendd(:,:,1,:)=U(:,0:nxg1,     1,0:nzg1)
    sendb(:,:,:,1)=U(:,0:nxg1,0:nyg1,    nz)
    sendf(:,:,:,1)=U(:,0:nxg1,0:nyg1,     1)

!    print*, 'Arrays allocated in rank ', rank
    
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
    !! Boundaries in yz faces
    !! Left boundary
    if (left>=0) then
       U(:,0,0:nyg1,0:nzg1) = recvl(:,1,:,:)
    else
       if (LB == BC_OUTFLOW) then
          U(:  ,0,0:nyg1,0:nzg1) =  U(:  ,1,0:nyg1,0:nzg1)
       elseif (LB == BC_REFLECTIVE) then
          U(1  ,0,0:nyg1,0:nzg1) =  U(1  ,1,0:nyg1,0:nzg1)
          U(2  ,0,0:nyg1,0:nzg1) = -U(2  ,1,0:nyg1,0:nzg1)
          U(3:5,0,0:nyg1,0:nzg1) =  U(3:5,1,0:nyg1,0:nzg1)
       elseif (LB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
!    print*, 'Left boundaries set in rank ', rank
    !! Right boundary
    if (right>=0) then
       U(:,nxg1,0:nyg1,0:nzg1) = recvr(:,1,:,:)
    else
       if (RB == BC_OUTFLOW) then
          U(:  ,nxg1,0:nyg1,0:nzg1) =  U(:  ,nx,0:nyg1,0:nzg1)
       elseif (RB == BC_REFLECTIVE) then
          U(1  ,nxg1,0:nyg1,0:nzg1) =  U(1  ,nx,0:nyg1,0:nzg1)
          U(2  ,nxg1,0:nyg1,0:nzg1) = -U(2  ,nx,0:nyg1,0:nzg1)
          U(3:5,nxg1,0:nyg1,0:nzg1) =  U(3:5,nx,0:nyg1,0:nzg1)
       elseif (RB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
!    print*, 'Right boundaries set in rank ', rank
    !! Boundaries in xz faces
    !! Down boundary
    if (down>=0) then
       U(:,0:nxg1,0,0:nzg1) = recvd(:,:,1,:)
    else
       if (DB == BC_OUTFLOW) then
          U(:  ,0:nxg1,0,0:nzg1) =  U(:  ,0:nxg1,1,0:nzg1)
       elseif (DB == BC_REFLECTIVE) then
          U(1:2,0:nxg1,0,0:nzg1) =  U(1:2,0:nxg1,1,0:nzg1)
          U(3  ,0:nxg1,0,0:nzg1) = -U(3  ,0:nxg1,1,0:nzg1)
          U(4:5,0:nxg1,0,0:nzg1) =  U(4:5,0:nxg1,1,0:nzg1)
       elseif (DB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
!    print*, 'Down boundaries set in rank ', rank
    !! Top boundary
    if (top>=0) then
       U(:,0:nxg1,nyg1,0:nzg1) = recvt(:,:,1,:)
    else
       if (TB == BC_OUTFLOW) then
          U(:  ,0:nxg1,nyg1,0:nzg1) =  U(:  ,0:nxg1,ny,0:nzg1)
       elseif (TB == BC_REFLECTIVE) then
          U(1:2,0:nxg1,nyg1,0:nzg1) =  U(1:2,0:nxg1,ny,0:nzg1)
          U(3  ,0:nxg1,nyg1,0:nzg1) = -U(3  ,0:nxg1,ny,0:nzg1)
          U(4:5,0:nxg1,nyg1,0:nzg1) =  U(4:5,0:nxg1,ny,0:nzg1)
       elseif (TB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif   
!    print*, 'Top boundaries set in rank ', rank
    !! Boundaries in xy faces
    !! Front boundary
    if (front>=0) then
       U(:,0:nxg1,0:nyg1,0) = recvf(:,:,:,1)
    else
       if (FB == BC_OUTFLOW) then
          U(:  ,0:nxg1,0:nyg1,0) =  U(:  ,0:nxg1,0:nyg1,1)
       elseif (FB == BC_REFLECTIVE) then
          U(1:3,0:nxg1,0:nyg1,0) =  U(1:3,0:nxg1,0:nyg1,1)
          U(4  ,0:nxg1,0:nyg1,0) = -U(4  ,0:nxg1,0:nyg1,1)
          U(5  ,0:nxg1,0:nyg1,0) =  U(5  ,0:nxg1,0:nyg1,1)
       elseif (FB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
!    print*, 'Front boundaries set in rank ', rank
    !! Back boundary
    if (back>=0) then
       U(:,0:nxg1,0:nyg1,nzg1) = recvb(:,:,:,1)
    else
       if (BB == BC_OUTFLOW) then
          U(:  ,0:nxg1,0:nyg1,nzg1) =  U(:  ,0:nxg1,0:nyg1,nz)
       elseif (BB == BC_REFLECTIVE) then
          U(1:3,0:nxg1,0:nyg1,nzg1) =  U(1:3,0:nxg1,0:nyg1,nz)
          U(4  ,0:nxg1,0:nyg1,nzg1) = -U(4  ,0:nxg1,0:nyg1,nz)
          U(5  ,0:nxg1,0:nyg1,nzg1) =  U(5  ,0:nxg1,0:nyg1,nz)
       elseif (BB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
!    print*, 'Back boundaries set in rank ', rank
#else
    
    !! Boundaries in yz faces
    !! Left boundary
    if (LB == BC_OUTFLOW) then
       U(:  ,0,:,:) =  U(:  ,1 ,:,:)
    elseif (LB == BC_REFLECTIVE) then
       U(1  ,0,:,:) =  U(1  ,1 ,:,:)
       U(2  ,0,:,:) = -U(2  ,1 ,:,:)
       U(3:5,0,:,:) =  U(3:5,1 ,:,:)
    elseif (LB == BC_PERIODIC) then
       U(:  ,0,:,:) =  U(:  ,nx,:,:)
    elseif (LB == BC_INFLOW) then
       !U(:,0,:,:) = ULin(:,:,:,:)
    endif
    !! Right boundary
    if (RB == BC_OUTFLOW) then
       U(:  ,nxg1,:,:) =  U(:  ,nx,:,:)
    elseif (RB == BC_REFLECTIVE) then
       U(1  ,nxg1,:,:) =  U(1  ,nx,:,:)
       U(2  ,nxg1,:,:) = -U(2  ,nx,:,:)
       U(3:5,nxg1,:,:) =  U(3:5,nx,:,:)
    elseif (RB == BC_PERIODIC) then
       U(:  ,nxg1,:,:) =  U(:  ,1 ,:,:)
    elseif (RB == BC_INFLOW) then
       !U(:,nxg1,:,:) = URin(:,:,:,:)
    endif

    !! Boundaries in xz faces
    !! Top boundary
    if (TB == BC_OUTFLOW) then
       U(:  ,:,nyg1,:) =  U(:  ,:,ny,:)
    elseif (TB == BC_REFLECTIVE) then
       U(1:2,:,nyg1,:) =  U(1:2,:,ny,:)
       U(3  ,:,nyg1,:) = -U(3  ,:,ny,:)
       U(4:5,:,nyg1,:) =  U(4:5,:,ny,:)
    elseif (TB == BC_PERIODIC) then
       U(:  ,:,nyg1,:) =  U(:  ,:,1 ,:)
    elseif (TB == BC_INFLOW) then
       !U(:,:,nyg1,:) = UTin(:,:,:,:)
    endif
    !! Down boundary
    if (DB == BC_OUTFLOW) then
       U(:  ,:,0,:) =  U(:  ,:,1 ,:)
    elseif (DB == BC_REFLECTIVE) then
       U(1:2,:,0,:) =  U(1:2,:,1 ,:)
       U(3  ,:,0,:) = -U(3  ,:,1 ,:)
       U(4:5,:,0,:) =  U(4:5,:,1 ,:)
    elseif (DB == BC_PERIODIC) then
       U(:  ,:,0,:) =  U(:  ,:,ny,:)
    elseif (DB == BC_INFLOW) then
       !U(:,:,0,:) = UDin(:,:,:,:)
    endif
        
    !! Boundaries in xy faces
    !! Back boundary
    if (BB == BC_OUTFLOW) then
       U(:  ,:,:,nzg1) =  U(:  ,:,:,nz)
    elseif (BB == BC_REFLECTIVE) then
       U(1:3,:,:,nzg1) =  U(1:3,:,:,nz)
       U(4  ,:,:,nzg1) = -U(4  ,:,:,nz)
       U(5  ,:,:,nzg1) =  U(5  ,:,:,nz)
    elseif (BB == BC_PERIODIC) then
       U(:  ,:,:,nzg1) =  U(:  ,:,:,1 )
    elseif (BB == BC_INFLOW) then
       !U(:,:,:,nzg1) = UBin(:,:,:,:)
    endif
    !! Front boundary
    if (FB == BC_OUTFLOW) then
       U(:  ,:,:,0) =  U(:  ,:,:,1 )
    elseif (FB == BC_REFLECTIVE) then
       U(1:3,:,:,0) =  U(1:3,:,:,1 )
       U(4  ,:,:,0) = -U(4  ,:,:,1 )
       U(5  ,:,:,0) =  U(5  ,:,:,1 )
    elseif (FB == BC_PERIODIC) then
       U(:  ,:,:,0) =  U(:  ,:,:,nz)
    elseif (FB == BC_INFLOW) then
       !U(:,:,:,0) = UFin(:,:,:,:)
    endif

#endif
    
  end subroutine boundariesI
  
  subroutine boundariesII(U)
    
    use parameters, only: nx, ny, nz, nxmin, nymin, nzmin, nxmax, nymax, nzmax, nghost, neq, np, &
         LB, RB, TB, DB, BB, FB
#ifdef MPIP
    use parameters, only: mpi_real_kind
    use globals, only: rank, err, left, right, down, top, front, back, comm3d
#endif
    use constants
    
    implicit none
    
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(inout) :: U
    
#ifdef MPIP
    
    include "mpif.h"
    
    integer                                        ::  status(MPI_STATUS_SIZE)
    real, dimension(neq,2,nymin:nymax,nzmin:nzmax) ::  sendr,recvr,sendl,recvl
    real, dimension(neq,nxmin:nxmax,2,nzmin:nzmax) ::  sendt,recvt,sendd,recvd
    real, dimension(neq,nxmin:nxmax,nymin:nymax,2) ::  sendb,recvb,sendf,recvf
    integer, parameter                             ::  srsizex = neq*nghost*(ny+2*nghost)*(nz+2*nghost)
    integer, parameter                             ::  srsizey = neq*nghost*(nx+2*nghost)*(nz+2*nghost)
    integer, parameter                             ::  srsizez = neq*nghost*(nx+2*nghost)*(ny+2*nghost)
#endif
    integer, parameter                             ::  nxm1 = nx-1, nym1 = ny-1, nzm1 = nz-1
    
#ifdef MPIP
    
    sendr(:,1:2,:,:)=U(:,nxm1:nx,:,:)
    sendl(:,1:2,:,:)=U(:,1:2    ,:,:)
    sendt(:,:,1:2,:)=U(:,:,nym1:ny,:)
    sendd(:,:,1:2,:)=U(:,:,1:2    ,:)
    sendb(:,:,:,1:2)=U(:,:,:,nzm1:nz)
    sendf(:,:,:,1:2)=U(:,:,:,1:2    )

    call mpi_sendrecv(sendr, srsizex, mpi_real_kind, right, 0,        &
         recvl, srsizex, mpi_real_kind, left , 0, comm3d, status, err)
       
    call mpi_sendrecv(sendl, srsizex, mpi_real_kind, left , 0,        &
         recvr, srsizex, mpi_real_kind, right, 0, comm3d, status, err)

    call mpi_sendrecv(sendt, srsizey, mpi_real_kind, top  , 0,        &
         recvd, srsizey, mpi_real_kind, down , 0, comm3d, status, err)
       
    call mpi_sendrecv(sendd, srsizey, mpi_real_kind, down , 0,        &
         recvt, srsizey, mpi_real_kind, top  , 0, comm3d, status, err)

    call mpi_sendrecv(sendb, srsizez, mpi_real_kind, back , 0,        &
         recvf, srsizez, mpi_real_kind, front, 0, comm3d, status, err)
       
    call mpi_sendrecv(sendf, srsizez, mpi_real_kind, front, 0,        &
         recvb, srsizez, mpi_real_kind, back , 0, comm3d, status, err)
    
    !! Boundaries in yz faces
    !! Left boundary
    if (left>=0) then
       U(:,-1:0,:,:) = recvl(:,1:2,:,:)
    else
       if (LB == BC_OUTFLOW) then
          U(:  ,0 ,:,:) =  U(:  ,1,:,:)
          U(:  ,-1,:,:) =  U(:  ,1,:,:)          
       elseif (LB == BC_REFLECTIVE) then
          U(1  ,0 ,:,:) =  U(1  ,1,:,:)
          U(2  ,0 ,:,:) = -U(2  ,1,:,:)
          U(3:5,0 ,:,:) =  U(3:5,1,:,:)
          U(1  ,-1,:,:) =  U(1  ,2,:,:)
          U(2  ,-1,:,:) = -U(2  ,2,:,:)
          U(3:5,-1,:,:) =  U(3:5,2,:,:)          
       elseif (LB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
    !! Right boundary
    if (right>=0) then
       U(:,nx+1:nx+2,:,:) = recvr(:,1:2,:,:)
    else
       if (RB == BC_OUTFLOW) then
          U(:  ,nx+1,:,:) =  U(:  ,nx,:,:)
          U(:  ,nx+2,:,:) =  U(:  ,nx,:,:)          
       elseif (RB == BC_REFLECTIVE) then
          U(1  ,nx+1,:,:) =  U(1  ,nx  ,:,:)
          U(2  ,nx+1,:,:) = -U(2  ,nx  ,:,:)
          U(3:5,nx+1,:,:) =  U(3:5,nx  ,:,:)
          U(1  ,nx+2,:,:) =  U(1  ,nxm1,:,:)
          U(2  ,nx+2,:,:) = -U(2  ,nxm1,:,:)
          U(3:5,nx+2,:,:) =  U(3:5,nxm1,:,:)
       elseif (RB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif

    !! Boundaries in xz faces
    !! Down boundary
    if (down>=0) then
       U(:,:,-1:0,:) = recvd(:,:,1:2,:)
    else
       if (DB == BC_OUTFLOW) then
          U(:  ,:, 0,:) =  U(:  ,:,1,:)
          U(:  ,:,-1,:) =  U(:  ,:,1,:)          
       elseif (DB == BC_REFLECTIVE) then
          U(1:2,:, 0,:) =  U(1:2,:,1,:)
          U(3  ,:, 0,:) = -U(3  ,:,1,:)
          U(4:5,:, 0,:) =  U(4:5,:,1,:)
          U(1:2,:,-1,:) =  U(1:2,:,2,:)
          U(3  ,:,-1,:) = -U(3  ,:,2,:)
          U(4:5,:,-1,:) =  U(4:5,:,2,:)          
       elseif (DB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
    !! Top boundary
    if (top>=0) then
       U(:,:,ny+1:ny+2,:) = recvt(:,:,1:2,:)
    else
       if (TB == BC_OUTFLOW) then
          U(:  ,:,ny+1,:) =  U(:  ,:,ny,:)
          U(:  ,:,ny+2,:) =  U(:  ,:,ny,:)          
       elseif (TB == BC_REFLECTIVE) then
          U(1:2,:,ny+1,:) =  U(1:2,:,ny  ,:)
          U(3  ,:,ny+1,:) = -U(3  ,:,ny  ,:)
          U(4:5,:,ny+1,:) =  U(4:5,:,ny  ,:)
          U(1:2,:,ny+2,:) =  U(1:2,:,nym1,:)
          U(3  ,:,ny+2,:) = -U(3  ,:,nym1,:)
          U(4:5,:,ny+2,:) =  U(4:5,:,nym1,:)          
       elseif (TB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif   

    !! Boundaries in xy faces
    !! Front boundary
    if (front>=0) then
       U(:,:,:,-1:0) = recvf(:,:,:,1:2)
    else
       if (FB == BC_OUTFLOW) then
          U(:  ,:,:, 0) =  U(:  ,:,:,1)
          U(:  ,:,:,-1) =  U(:  ,:,:,1)          
       elseif (FB == BC_REFLECTIVE) then
          U(1:3,:,:, 0) =  U(1:3,:,:,1)
          U(4  ,:,:, 0) = -U(4  ,:,:,1)
          U(5  ,:,:, 0) =  U(5  ,:,:,1)
          U(1:3,:,:,-1) =  U(1:3,:,:,2)
          U(4  ,:,:,-1) = -U(4  ,:,:,2)
          U(5  ,:,:,-1) =  U(5  ,:,:,2)
       elseif (FB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif
    !! Back boundary
    if (back>=0) then
       U(:,:,:,nz+1:nz+2) = recvb(:,:,:,1:2)
    else
       if (BB == BC_OUTFLOW) then
          U(:  ,:,:,nz+1) =  U(:  ,:,:,nz)
          U(:  ,:,:,nz+2) =  U(:  ,:,:,nz)          
       elseif (BB == BC_REFLECTIVE) then
          U(1:3,:,:,nz+1) =  U(1:3,:,:,nz  )
          U(4  ,:,:,nz+1) = -U(4  ,:,:,nz  )
          U(5  ,:,:,nz+1) =  U(5  ,:,:,nz  )
          U(1:3,:,:,nz+2) =  U(1:3,:,:,nzm1)
          U(4  ,:,:,nz+2) = -U(4  ,:,:,nzm1)
          U(5  ,:,:,nz+2) =  U(5  ,:,:,nzm1)
       elseif (BB == BC_INFLOW) then
          !U(:,0,0:nyg1,0:nzg1) = ULin(:,0,:,:)
       endif
    endif

#else
    !! boundaries in zy faces
    if (LB == BC_PERIODIC) then
       U(:,nxmin:0,:,:) = U(:,nxm1:nx,:,:)
    elseif (LB == BC_OUTFLOW) then
       U(:,0,:,:) = U(:,1,:,:)
       U(:,nxmin,:,:) = U(:,1,:,:)
    elseif (LB == BC_REFLECTIVE) then
       U(1,0,:,:) = U(1,1,:,:)
       U(2,0,:,:) = -U(2,1,:,:)
       U(3:5,0,:,:) = U(3:5,1,:,:)
       U(1,nxmin,:,:) = U(1,nghost,:,:)
       U(2,nxmin,:,:) = -U(2,nghost,:,:)
       U(3:5,nxmin,:,:) = U(3:5,nghost,:,:)
    elseif (LB == BC_INFLOW) then
       !U(:,nxmin:0,:,:) = ULin(:,1:2,:,:)
    endif
    
    if (RB == BC_PERIODIC) then
       U(:,nx+1:nxmax,:,:) = U(:,1:nghost,:,:)
    elseif (RB == BC_OUTFLOW) then
       U(:,nx+1,:,:) = U(:,nx,:,:)
       U(:,nxmax,:,:) = U(:,nx,:,:)
    elseif (RB == BC_REFLECTIVE) then
       U(1,nx+1,:,:) = U(1,nx,:,:)
       U(2,nx+1,:,:) = -U(2,nx,:,:)
       U(3:5,nx+1,:,:) = U(3:5,nx,:,:)
       U(1,nxmax,:,:) = U(1,nxm1,:,:)
       U(2,nxmax,:,:) = -U(2,nxm1,:,:)
       U(3:5,nxmax,:,:) = U(3:5,nxm1,:,:)
    elseif (RB == BC_INFLOW) then
       !U(:,nx+1:nxmax,:,:) = URin(:,:,:,:)
    endif
        
    !! boundaries in xz faces
    
    if (DB == BC_REFLECTIVE) then
       U(1:2,:,0,:) = U(1:2,:,1,:)
       U(1:2,:,-1,:) = U(1:2,:,2,:)
       U(3,:,0,:) = -U(3,:,1,:)
       U(3,:,-1,:) = -U(3,:,2,:)
       U(4:5,:,0,:) = U(4:5,:,1,:)
       U(4:5,:,-1,:) = U(4:5,:,2,:)
    elseif (DB == BC_OUTFLOW) then
       U(:,:,0,:) = U(:,:,1,:)
       U(:,:,-1,:) = U(:,:,1,:)
    elseif (DB == BC_PERIODIC) then
       U(:,:,nymin:0,:) = U(:,:,ny-1:ny,:)
    elseif (DB == BC_INFLOW) then
       !U(:,:,nymin:0,:) = UDin(:,:,:,:)
    endif
    
    if (TB == BC_REFLECTIVE) then
       U(1:2,:,ny+1,:) = U(1:2,:,ny,:)
       U(1:2,:,ny+2,:) = U(1:2,:,ny-1,:)
       U(3,:,ny+1,:) = -U(3,:,ny,:)
       U(3,:,ny+2,:) = -U(3,:,ny-1,:)
       U(4:5,:,ny+1,:) = U(4:5,:,ny,:)
       U(4:5,:,ny+2,:) = U(4:5,:,ny-1,:)
    elseif (TB == BC_OUTFLOW) then
       U(:,:,ny+1,:) = U(:,:,ny,:)
       U(:,:,ny+2,:) = U(:,:,ny,:)
    elseif (TB == BC_PERIODIC) then
       U(:,:,ny+1:nymax,:) = U(:,:,1:2,:)
    elseif (TB == BC_INFLOW) then
       !U(:,:,ny+1:nymax,:) = UTin(:,:,:,:)
    endif

    !! boundaries in xy faces
    if (FB == BC_REFLECTIVE) then
       U(1:3,:,:,0) = U(1:3,:,:,1)
       U(1:3,:,:,-1) = U(1:3,:,:,2)
       U(4,:,:,0) = -U(4,:,:,1)
       U(4,:,:,-1) = -U(4,:,:,2)
       U(5,:,:,0) = U(5,:,:,1)
       U(5,:,:,-1) = U(5,:,:,2)
    elseif (FB == BC_OUTFLOW) then
       U(:,:,:,0) = U(:,:,:,1)
       U(:,:,:,-1) = U(:,:,:,1)
    elseif (FB == BC_PERIODIC) then
       U(:,:,:,nzmin:0) = U(:,:,:,nz-1:nz)
    elseif (FB == BC_INFLOW) then
       !U(:,:,:,nzmin:0) = UFin(:,:,:,:)
    endif
    
    if (BB == BC_REFLECTIVE) then
       U(1:3,:,:,nz+1) = U(1:3,:,:,nz)
       U(1:3,:,:,nz+2) = U(1:3,:,:,nz-1)
       U(4,:,:,nz+1) = -U(4,:,:,nz)
       U(4,:,:,nz+2) = -U(4,:,:,nz-1)
       U(5,:,:,nz+1) = U(5,:,:,nz)
       U(5,:,:,nz+2) = U(5,:,:,nz-1)
    elseif (BB == BC_OUTFLOW) then
       U(:,:,:,nz+1) = U(:,:,:,nz)
       U(:,:,:,nz+2) = U(:,:,:,nz)
    elseif (BB == BC_PERIODIC) then
       U(:,:,:,nz+1:nzmax) = U(:,:,:,1:2)
    elseif (BB == BC_INFLOW) then
       !U(:,:,:,nz+1:nzmax) = UBin(:,:,:,:)
    endif

#endif
    
  end subroutine boundariesII
