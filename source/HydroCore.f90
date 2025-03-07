  subroutine calc_prims(UPPrim,UU)

    use parameters, only: nxmin, nymin, nzmin, nxmax, nymax, nzmax, gamma, neq

    implicit none

    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(in) :: UU
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(out) :: UPPrim
    integer                         :: i,j,k!,s
    
    do i = nxmin, nxmax
       do j = nymin, nymax
          do k = nzmin, nzmax
             UPPrim(1,i,j,k) = UU(1,i,j,k)
             UPPrim(2,i,j,k) = UU(2,i,j,k)/UU(1,i,j,k)
             UPPrim(3,i,j,k) = UU(3,i,j,k)/UU(1,i,j,k)
             UPPrim(4,i,j,k) = UU(4,i,j,k)/UU(1,i,j,k)
             UPPrim(5,i,j,k) = (gamma-1.0)*(UU(5,i,j,k) - 0.5*(UU(2,i,j,k)**2 + &
                  UU(3,i,j,k)**2 + UU(4,i,j,k)**2)/UU(1,i,j,k))

             UPPrim(5,i,j,k) = max(UPPrim(5,i,j,k),1e-20)
!             if (UPPrim(5,i,j,k).ne.UPPrim(5,i,j,k)) then
!                s = s + 1
!                UPPrim(5,i,j,k) = 0.001
!                print*, 'Pressure Neg/NaNs'
!                stop
!             endif
          enddo
       enddo
    enddo
!    print*, 's = ', s
    
  end subroutine calc_prims
  
  subroutine xcoord(i, x)

    use globals, only: coords
    use parameters, only: xl, xr, mpix, dx
    
    integer, intent(in)     ::  i
    real, intent(out)       ::  x
    
    x = xl + coords(0)*(xr-xl)/mpix + dx*i - dx/2.0
    
    return
    
  end subroutine xcoord
  
  subroutine ycoord(j, y)

    use globals, only: coords
    use parameters, only: yl, yr, dy, mpiy
    
    integer, intent(in)     ::  j
    real, intent(out)       ::  y
    
    y = yl + coords(1)*(yr-yl)/mpiy + dy*j - dy/2.0
    
    return
    
  end subroutine ycoord
  
  subroutine zcoord(k, z)

    use globals, only: coords
    use parameters, only: zl, zr, dz, mpiz
    
    integer, intent(in)     ::  k
    real, intent(out)       ::  z
    
    z = zl + coords(2)*(zr-zl)/mpiz + dz*k - dz/2.0
    
    return
    
  end subroutine zcoord

  subroutine iindex(x, i)

    use globals, only: coords
    use parameters, only: xl, xr, dx, mpix

    integer, intent(out)    :: i
    real, intent(in)        :: x

    i = nint((x - xl - coords(0)*(xr-xl)/mpix + dx/2.0)/dx)

  end subroutine iindex

  subroutine jindex(y, j)

    use globals, only: coords
    use parameters, only: yl, yr, dy, mpiy

    integer, intent(out)    :: j
    real, intent(in)        :: y

    j = nint((y - yl - coords(1)*(yr-yl)/mpiy + dy/2.0)/dy)

  end subroutine jindex

  subroutine kindex(z, k)

    use globals, only: coords
    use parameters, only: zl, zr, dz, mpiz

    integer, intent(out)    :: k
    real, intent(in)        :: z

    k = nint((z - zl - coords(2)*(zr-zl)/mpiz + dz/2.0)/dz)

  end subroutine kindex
  
  subroutine time_step(UPPrim, dt)
    
    use parameters, only: nx, ny, nz, dx, dy, dz, mpi_real_kind, cfl, nxmin, nxmax, &
         nymin, nymax, nzmin, nzmax, neq
    use globals, only: err, comm3d
    
    implicit none
#ifdef MPIP
    include "mpif.h"
#endif
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(in) :: UPPrim
    real, intent(out)       :: dt
    real                    :: dtp
    real                    :: maxvx, maxvy, maxvz, cs
    integer                 :: i,j,k
    
    maxvx = 0
    maxvy = 0
    maxvz = 0
    cs = 0
    
    do i = 1, nx
       do j = 1, ny
          do k = 1, nz
             call csound(UPPrim(5,i,j,k),UPPrim(1,i,j,k), cs)
             maxvx = max(maxvx, abs(UPPrim(2,i,j,k)) + cs)
             maxvy = max(maxvy, abs(UPPrim(3,i,j,k)) + cs)
             maxvz = max(maxvz, abs(UPPrim(4,i,j,k)) + cs)
          enddo
       enddo
    enddo
    dtp = cfl*min(dx/maxvx,dy/maxvy,dz/maxvz)/sqrt(3.0)
#ifdef MPIP
    call mpi_allreduce(dtp, dt, 1, mpi_real_kind, mpi_min, comm3d, err)
#else
    dt = dtp
#endif
    !if (rank == 0) then
    !write(*,*) 'dx = ', dx, ' maxvx = ', maxvx, ' dy = ', dy, ' maxvy = ', maxvy, ' dz = ', dz, ' maxvz = ', maxvz
    !write(*,*) 'dtp = ', dtp, 'dt = ', dt, 'cs = ', cs
    !write(*,*) maxval(UPPrim(2,:,:,:)), maxval(UPPrim(3,:,:,:)), maxval(UPPrim(4,:,:,:))
    !endif
  end subroutine time_step
  
  subroutine csound(p,den,cs)
    
    use parameters, only: gamma
    
    implicit none
    
    real, intent(in)     :: p, den
    real, intent(out)    :: cs
    
    cs = sqrt(gamma*p/den)
    
  end subroutine csound
  
  subroutine step(t,it,dt)
    
    use parameters, only: nx, ny, nz, neq, eta, nxmax, nymax, nzmax
    use globals, only: U, UP
    
    implicit none
    
    integer                     :: i,j,k,e
    integer, intent(inout)      :: it
    real, intent(inout)         :: t
    real, intent(in)            :: dt
    
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             do e = 1,neq
                U(e,i,j,k) = UP(e,i,j,k) + eta*(U(e,i+1,j,k) + U(e,i-1,j,k) + U(e,i,j+1,k) + &
                     U(e,i,j-1,k) + U(e,i,j,k+1) + U(e,i,j,k-1) - 6.0*U(e,i,j,k))
             enddo
          enddo
       enddo
    enddo
    
    U(:,0,:,:) = UP(:,0,:,:)
    U(:,-1,:,:) = UP(:,-1,:,:)
    U(:,nxmax-1,:,:) = UP(:,nxmax-1,:,:)
    U(:,nxmax,:,:) = UP(:,nxmax,:,:)
    U(:,:,0,:) = UP(:,:,0,:)
    U(:,:,-1,:) = UP(:,:,-1,:)
    U(:,:,nymax-1,:) = UP(:,:,nymax-1,:)
    U(:,:,nymax,:) = UP(:,:,nymax,:)
    U(:,:,:,0) = UP(:,:,:,0)
    U(:,:,:,-1) = UP(:,:,:,-1)
    U(:,:,:,nzmax-1) = UP(:,:,:,nzmax-1)
    U(:,:,:,nzmax) = UP(:,:,:,nzmax)
    
    t = t + dt
    it = it + 1
    
  end subroutine step
  
  subroutine prim2f(ff,prim,dim)
    
    use parameters, only: gamma, neq
    
    implicit none
    
    real, dimension(neq), intent(out)         :: ff
    real, dimension(neq), intent(in)          :: prim
    integer, intent(in)                       :: dim
    
    select case(dim)
       
    case(1)
       
       ff(1) = prim(1)*prim(2)
       ff(2) = prim(1)*prim(2)*prim(2) + prim(5)
       ff(3) = prim(1)*prim(2)*prim(3)
       ff(4) = prim(1)*prim(2)*prim(4)
       ff(5) = prim(2)*(0.5*prim(1)*(prim(2)**2 + prim(3)**2 + prim(4)**2) + (gamma/(gamma-1))*prim(5))
       
    case(2)
       
       ff(1) = prim(1)*prim(3)
       ff(2) = prim(1)*prim(2)*prim(3)
       ff(3) = prim(1)*prim(3)*prim(3) + prim(5)
       ff(4) = prim(1)*prim(3)*prim(4)
       ff(5) = prim(3)*(0.5*prim(1)*(prim(2)**2 + prim(3)**2 + prim(4)**2) + (gamma/(gamma-1))*prim(5))
       
    case(3)
       
       ff(1) = prim(1)*prim(4)
       ff(2) = prim(1)*prim(2)*prim(4)
       ff(3) = prim(1)*prim(3)*prim(4)
       ff(4) = prim(1)*prim(4)*prim(4) + prim(5)
       ff(5) = prim(4)*(0.5*prim(1)*(prim(2)**2 + prim(3)**2 + prim(4)**2) + (gamma/(gamma-1))*prim(5))
       
    end select
    
  end subroutine prim2f
  
  subroutine output(nout, tout)
    
    use parameters, only: nxmin, nymin, nzmin, nxmax, nymax, nzmax, nghost, dx, dy, dz, dtout, outfile, &
         xr, xl, yr, yl, zr, zl, np, mpix, mpiy, mpiz, outputpath
    use globals, only: UPrim, rank, coords
#ifdef GRAV
    use globals, only: PHI
#endif
    use constants, only: DAT, VTK, BIN, PC
    implicit none
    
    integer, intent(inout) :: nout
    real, intent(inout)    :: tout
    integer                :: i, j, k, unitout
    character (len=128)    :: file1, file2, file3
    character(len=200)     :: cbuffer
    character(1)           :: lf = char(10)

    write(file3,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'points',rank,'N',nout,'.bin'    
    write(file1,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'points',rank,'N',nout,'.dat'
    write(file2,'(a,i3.3,a,i3.3,a)')  trim(outputpath)//'EPV',rank,'N',nout,'.vtk'
    unitout=10
    select case(outfile)

    case(BIN)

       open(unit=unitout,file=file3,status='unknown',access='stream',form='unformatted',convert='LITTLE_ENDIAN')

       write(unitout) UPrim(:,:,:,:)
       
       close(unitout)
#ifdef GRAV
       open(unit=unitout,file='phi'//file3,status='unknown',access='stream',form='unformatted',convert='LITTLE_ENDIAN')

       write(unitout) PHI(:,:,:)

       close(unitout)
#endif
    case(DAT)
       
       open(unit=unitout,file=file1,status='unknown')

#ifdef GRAV
       do k = nxmin+nghost, nxmax-nghost
          do j = nymin+nghost, nymax-nghost
             do i = nzmin+nghost, nzmax-nghost
                write(unitout,'(14(e16.7e3))') float(i)*dx,float(j)*dy,float(k)*dz,UPrim(:,i,j,k),phi(i,j,k)
             enddo
          enddo
       enddo
#else
       do k = nxmin+nghost, nxmax-nghost
          do j = nymin+nghost, nymax-nghost
             do i = nzmin+nghost, nzmax-nghost
                write(unitout,'(14(e16.7e3))') float(i)*dx,float(j)*dy,float(k)*dz,UPrim(:,i,j,k)
             enddo
          enddo
       enddo
#endif
       close(unitout)

    case(VTK)
       
       open(unit=unitout,file=file2,status='unknown',access='stream',form='unformatted',convert='BIG_ENDIAN')
       
       write(unitout) '# vtk DataFile Version 2.0', lf
       write(unitout) 'output from Diable', lf
       write(unitout) 'BINARY', lf
       write(unitout) 'DATASET STRUCTURED_POINTS', lf
       write(cbuffer, '("DIMENSIONS ",(3i6,1x))') nxmax-nghost,nymax-nghost,nzmax-nghost
       write(unitout) trim(cbuffer), lf
       
       write(cbuffer, '("ORIGIN "    ,3e15.7)') xl+coords(0)*(xr-xl)/mpix,yl+coords(1)*(yr-yl)/mpiy,zl+coords(2)*(zr-zl)/mpiz
       write(unitout) trim(cbuffer), lf
       !
       write(cbuffer, '("SPACING",3e15.7)') dx,dy,dz
       write(unitout) trim(cbuffer), lf
       !
       !   writes the variables, scalars first then vectors
       !
       write(cbuffer,'(a,i10)') 'POINT_DATA ', (nxmax-nghost)*(nymax-nghost)*(nzmax-nghost)
       write(unitout) trim(cbuffer), lf
       !
       write(unitout) 'SCALARS Density double 1', lf
       !
       write(unitout) 'LOOKUP_TABLE default', lf
       !
       do k = 1,nzmax-nghost
          do j = 1,nymax-nghost
             do i = 1,nxmax-nghost
                write(unitout) UPrim(1,i,j,k)
             enddo
          enddo
       enddo
       !
       write(unitout) lf
       write(unitout) 'SCALARS Pressure double 1', lf
       !
       write(unitout) 'LOOKUP_TABLE default', lf
       !
       do k = 1,nzmax-nghost
          do j = 1,nymax-nghost
             do i = 1,nxmax-nghost
                write(unitout) UPrim(5,i,j,k)
             enddo
          enddo
       enddo
#ifdef GRAV
       write(unitout) lf
       write(unitout) 'SCALARS Potential double 1', lf
       !
       write(unitout) 'LOOKUP_TABLE default', lf
       !
       do k = 1,nzmax-nghost
          do j = 1,nymax-nghost
             do i = 1,nxmax-nghost
                write(unitout) PHI(i,j,k)
             enddo
          enddo
       enddo
#endif
       write(unitout) lf
       write(unitout) 'VECTORS Velocity double', lf
       !
       do k = 1,nzmax-nghost
          do j = 1,nymax-nghost
             do i = 1,nxmax-nghost
                write(unitout) UPrim(2,i,j,k), UPrim(3,i,j,k), UPrim(4,i,j,k)
             enddo
          enddo
       enddo
       write(unitout) lf
       close(unitout)
       
    endselect

    if (rank == 0) then
       write(*,*) 'File ', nout, ' writen'
    endif
    nout = nout + 1
    tout = dtout*nout

  end subroutine output
  
