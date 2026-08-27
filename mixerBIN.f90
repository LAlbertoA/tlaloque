program mixer

  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONSTANTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                             !! Distance Scales
  real, parameter                            :: PC = 3.085677588e+18, AU = 1.495978707e+13, KM = 1.0e5, NO = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SIMULATION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                             !! GRAVITY = .true./.false. 
  logical, parameter                         :: GRAVITY = .true.
                                             !! mpix/y/z = [Number of cores in x/y/z], nx/y/z = [Number of Cells in x/y/z]
  integer, parameter                         :: mpix = 2, mpiy = 2, mpiz = 2, nx = 64, ny = 64, nz = 64
                                             !! nout = [Number of outputs per core]
  integer, parameter                         :: nout = 11
                                             !! x/y/zl = [minimum x/y/z value in domain], x/y/zr = [maximum x/y/z value in domain]
  real, parameter                            :: xl = -2, yl = -2, zl = -2, xr = 2, yr = 2, zr = 2
                                             !! scale = [PC/AU/KM/N Scale factor for the distances (cell size & origin)]
  real, parameter                            :: scale = NO
                                             !! folder = [Files location]
  character(*), parameter                    :: folder = '/home/luis-arcos/Documents/tlaloque/DATA/'
                                             !! Multiproc filenames (filenames should be in the form of prefname*nproc*'N'*nout*'.vtk')
  character(*), parameter                    :: prefname = 'EPV'
                                             !! Output filename (filename will be in the form of postfname*nout*'.vtk')
  character(*), parameter                    :: postfname = 'EVR'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DON'T MOVE ANYTHING BELOW THIS LINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter                         :: nprocs = mpix*mpiy*mpiz, unitout = 10
  integer                                    :: p, o, i, j, k, n
  integer, dimension(3,nprocs)               :: test
  real, dimension(3)                         :: origin, dxdydz
  real, dimension(:,:,:,:), allocatable      :: U
  real, dimension(:,:,:), allocatable        :: PHI

  character(len=230), dimension(nprocs,nout) :: filenames
  character(len=200)                         :: fname, fout, cbuffer
  character(1)                               :: lf = char(10)
  character                                  :: rbuf(283), rbuf2(48), rbuf3(49), rbuf4(25)

  dxdydz(1) = (xr-xl)/(nx*mpix)/scale
  dxdydz(2) = (yr-yl)/(ny*mpiy)/scale
  dxdydz(3) = (zr-zl)/(nz*mpiz)/scale

  origin(1) = xl/scale
  origin(2) = yl/scale
  origin(3) = zl/scale
  
  do p = 1,nprocs
     do o = 1,nout
        write(fname,'(a,i3.3,a,i3.3,a)') prefname//'',p-1,'N',o-1,'.vtk'
        filenames(p,o) = folder//fname
!        print*, filenames(p,o)
     enddo
  enddo

  do i = 1,mpix
     do j = 1,mpiy
        do k = 1,mpiz
           test(1,(i-1)*mpiz*mpiy+(j-1)*mpiz+k) = i-1
           test(2,(i-1)*mpiz*mpiy+(j-1)*mpiz+k) = j-1
           test(3,(i-1)*mpiz*mpiy+(j-1)*mpiz+k) = k-1
        enddo
     enddo
  enddo

  allocate(U(5,nx*mpix,ny*mpiy,nz*mpiz))
  allocate(PHI(nx*mpix,ny*mpiy,nz*mpiz))
  
  do n = 1,nout
     do p = 1,nprocs
        open(10,file=filenames(p,n),status='old',access='stream',form='unformatted',convert='BIG_ENDIAN')
        
        read(10) rbuf
        
        read(10) U(1,1+nx*test(1,p):nx+nx*test(1,p),1+ny*test(2,p):ny+ny*test(2,p),1+nz*test(3,p):nz+nz*test(3,p))
        
        read(10) rbuf2

        read(10) U(5,1+nx*test(1,p):nx+nx*test(1,p),1+ny*test(2,p):ny+ny*test(2,p),1+nz*test(3,p):nz+nz*test(3,p))

        if (GRAVITY .eqv. .true.) then
           read(10) rbuf3

           read(10) PHI(1+nx*test(1,p):nx+nx*test(1,p),1+ny*test(2,p):ny+ny*test(2,p),1+nz*test(3,p):nz+nz*test(3,p))
        endif
        read(10) rbuf4

        read(10) U(2:4,1+nx*test(1,p):nx+nx*test(1,p),1+ny*test(2,p):ny+ny*test(2,p),1+nz*test(3,p):nz+nz*test(3,p))

        close(10)
        print*, trim(filenames(p,n))//' file read'

     enddo
     
     write(fout,'(a,i3.3,a)') folder//postfname//'',n-1,'BIN.vtk'

     print'(a,a)', 'Writing ',fout
     open(unit=unitout,file=fout,status='new',access='stream',form='unformatted',convert='BIG_ENDIAN')
     
     write(unitout) '# vtk DataFile Version 2.0', lf
     write(unitout) 'output from Diable', lf
     write(unitout) 'BINARY', lf
     write(unitout) 'DATASET STRUCTURED_POINTS', lf
     write(cbuffer, '("DIMENSIONS ",(3i6,1x))') nx*mpix,ny*mpiy,nz*mpiz
     write(unitout) trim(cbuffer), lf
     
     write(cbuffer, '("ORIGIN "    ,3e15.7)') origin(1),origin(2),origin(3)
     write(unitout) trim(cbuffer), lf
     !
     
     write(cbuffer, '("SPACING",3e15.7)') dxdydz(1),dxdydz(2),dxdydz(3)
     write(unitout) trim(cbuffer), lf
     !                                                                                                                             
     !   writes the variables, scalars first then vectors                                                                          
     !                                                                                                                             
     write(cbuffer,'(a,i10)') 'POINT_DATA ', (nx*mpix)*(ny*mpiy)*(nz*mpiz)
     write(unitout) trim(cbuffer), lf
     !                                                                                                                             
     write(unitout) 'SCALARS Density double 1', lf
     !                                                                                                                             
     write(unitout) 'LOOKUP_TABLE default', lf
     !                                                                                                                             
     write(unitout) U(1,:,:,:)
     !
     write(unitout) lf
     write(unitout) 'SCALARS Pressure double 1', lf
     !                                                                                                                             
     write(unitout) 'LOOKUP_TABLE default', lf
     !                                                                                                                             
     write(unitout) U(5,:,:,:)
     !
     if (GRAVITY .eqv. .true.) then
        write(unitout) lf
        write(unitout) 'SCALARS Potential double 1', lf
        !                                                                                                                             
        write(unitout) 'LOOKUP_TABLE default', lf
        !                                                                                                                             
        write(unitout) PHI(:,:,:)
        !
     endif
     write(unitout) lf
     write(unitout) 'VECTORS Velocity double', lf
     !                                                                                                                             
     write(unitout) U(2:4,:,:,:)
     write(unitout) lf
     close(unitout)
  enddo  

end program mixer
