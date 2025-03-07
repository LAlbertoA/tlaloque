module HLL

contains

  subroutine Godunov(UPPP,UU,FF,GG,HH,dt)

    use parameters, only: dx, dy, dz, nx, ny, nz, neq, nxmin, nymin, nzmin, nxmax, nymax, nzmax
#ifdef GRAV
    use sources
    use globals, only: PHI, UPrim, S
#endif
    implicit none

    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(out) :: UPPP
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(in) :: UU, FF, GG, HH
    real, intent(in)                :: dt
    real, dimension(nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1) :: dens
    real                            :: dtdx, dtdy, dtdz
    integer                         :: i,j,k,e!,png
#ifdef MPIP
    integer                         :: err
#endif

    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz
#ifdef GRAV
    dens = UU(1,nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1)
    !dens = UU(1,nxmin+2:nxmax-2, nymin+2:nymax-2, nzmin+2:nzmax-2)

    !call poisson_solver1(dens, PHI, 3)
    call MultiGrid(dens,PHI)
    call S_Solver(UPrim,PHI,S)
#endif
!        print*, maxval(UPrim(2,:,:,:)),maxval(UPrim(3,:,:,:)),maxval(UPrim(4,:,:,:)),maxval(UPrim(5,:,:,:))
!        print*, 'Minimos UPrim', minval(UPrim(2,:,:,:)),minval(UPrim(3,:,:,:)),minval(UPrim(4,:,:,:)),minval(UPrim(5,:,:,:))
!        print*, dt
!    png = 0
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             do e = 1,neq

                UPPP(e,i,j,k) = UU(e,i,j,k) - dtdx*(FF(e,i,j,k) - FF(e,i-1,j,k)) - &
                     dtdy*(GG(e,i,j,k) - GG(e,i,j-1,k)) - dtdz*(HH(e,i,j,k) - HH(e,i,j,k-1))
#ifdef GRAV
                UPPP(e,i,j,k) = UPPP(e,i,j,k) + dt*S(e,i,j,k)
#endif

             !   if (UPPP(e,i,j,k).ne.UPPP(e,i,j,k)) then
             !      print*, 'Tienes NaNs', e,i,j,k,UPPP(e,i,j,k), UU(e,i,j,k), FF(e,i,j,k), GG(e,i,j,k),  HH(e,i,j,k)
             !      stop
             !   endif
             enddo
             !if (UPPP(5,i,j,k).ne.abs(UPPP(5,i,j,k))) then
             !   UPPP(5,i,j,k) = 1.0e-15
             !endif
!             if (UPPP(5,i,j,k).ne.abs(UPPP(5,i,j,k))) then
!                png = png + 1
                !print*, 'Pressure neg', UPPP(5,i,j,k), dtdx*(FF(1,i,j,k) - FF(1,i-1,j,k)), &
                !dtdy*(GG(1,i,j,k) - GG(1,i,j-1,k)), dtdz*(HH(1,i,j,k) - HH(1,i,j,k-1))
!             endif
          enddo
       enddo
    enddo
!    print*, 'Pressure negative/Nans', png
    !write(*,*) maxval(FF(:,:,:,:)), maxval(GG(:,:,:,:)), maxval(HH(:,:,:,:))
    
  end subroutine Godunov

  subroutine HLLFluxes(UU,FF,GG,HH,UPPrim,choice)

    use parameters, only: nxmin, nymin, nzmin, nx, ny, nz, neq, nxmax, nymax, nzmax, limtr
    use HydroCore, only: prim2f
    
    implicit none
    
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(in) :: UU, UPPrim
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(out) :: FF, GG, HH
    integer, intent(in)               :: choice
    real, dimension(neq)              :: primLL,primL,primR,primRR,fff, fl, fr, primU, primB
    real, dimension(neq)              :: primDD, primUU, primBB, primFF, primD, primF
    real                              :: sl, sr
    integer                           :: i,j,k!,clF,crF,cmF,clG,crG,cmG,clH,crH,cmH
!    clF = 0
!    crF = 0
!    cmF = 0
!    clG = 0
!    crG = 0
!    cmG = 0
!    clH = 0
!    crH = 0
!    cmH = 0
    
    !write(*,*) '!!!!GODUNOV DEBUG!!!!'
    !write(*,*) maxval(UPPrim)
    select case(choice)
       
    case(1)
       
       do i = nxmin+1, nx
          do j = nymin+1, ny
             do k = nzmin+1, nz
                primL(:) = UPPrim(:,i,j,k)
                primR(:) = UPPrim(:,i+1,j,k)
                call wavespeed(primL, primR, sl, sr, 1)
                if (sl>=0) then
                   call prim2f(fff, primL, 1)
                   FF(:,i,j,k) = fff(:)
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos L F'
                   !endif
                else if (sr<= 0) then
                   call prim2f(fff, primR, 1)
                   FF(:,i,j,k) = fff(:)
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos R F'
                   !endif
                else
                   call prim2f(fl, primL, 1)
                   call prim2f(fr, primR, 1)
                   FF(:,i,j,k) = (sr*fl(:)-sl*fr(:)+sl*sr*(UU(:,i+1,j,k)-UU(:,i,j,k)))/(sr-sl)
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos Mix F'
                   !endif
                endif
                primD(:) = UPPrim(:,i,j,k)
                primU(:) = UPPrim(:,i,j+1,k)
                call wavespeed(primD, primU, sl, sr, 2)
                if (sl>=0) then
                   call prim2f(fff, primD, 2)
                   GG(:,i,j,k) = fff(:)
                   !if (GG(1,i,j,k).ne.GG(1,i,j,k).or.GG(2,i,j,k).ne.GG(2,i,j,k).or.GG(3,i,j,k).ne.GG(3,i,j,k).or.&
                   !   GG(4,i,j,k).ne.GG(4,i,j,k).or.GG(5,i,j,k).ne.GG(5,i,j,k)) then
                   !   print*, 'Flujos L G'
                   !endif
                else if (sr<=0) then
                   call prim2f(fff, primU, 2)
                   GG(:,i,j,k) = fff(:)
                   !if (GG(1,i,j,k).ne.GG(1,i,j,k).or.GG(2,i,j,k).ne.GG(2,i,j,k).or.GG(3,i,j,k).ne.GG(3,i,j,k).or.&
                   !   GG(4,i,j,k).ne.GG(4,i,j,k).or.GG(5,i,j,k).ne.GG(5,i,j,k)) then
                   !   print*, 'Flujos R G'
                   !endif
                else
                   call prim2f(fl, primD, 2)
                   call prim2f(fr, primU, 2)
                   GG(:,i,j,k) = (sr*fl(:)-sl*fr(:)+sl*sr*(UU(:,i,j+1,k)-UU(:,i,j,k)))/(sr-sl)
                   !if (GG(1,i,j,k).ne.GG(1,i,j,k).or.GG(2,i,j,k).ne.GG(2,i,j,k).or.GG(3,i,j,k).ne.GG(3,i,j,k).or.&
                   !   GG(4,i,j,k).ne.GG(4,i,j,k).or.GG(5,i,j,k).ne.GG(5,i,j,k)) then
                   !   print*, 'Flujos Mix G'
                   !endif
                endif
                primF(:) = UPPrim(:,i,j,k)
                primB(:) = UPPrim(:,i,j,k+1)
                call wavespeed(primF, primB, sl, sr, 3)
                if (sl>=0) then
                   call prim2f(fff, primF, 3)
                   HH(:,i,j,k) = fff(:)
                   !if (HH(1,i,j,k).ne.HH(1,i,j,k).or.HH(2,i,j,k).ne.HH(2,i,j,k).or.HH(3,i,j,k).ne.HH(3,i,j,k).or.&
                   !   HH(3,i,j,k).ne.HH(3,i,j,k).or.HH(5,i,j,k).ne.HH(5,i,j,k)) then
                   !   print*, 'Flujos L H'
                   !endif
                else if (sr<=0) then
                   call prim2f(fff, primB, 3)
                   HH(:,i,j,k) = fff(:)
                   !if (HH(1,i,j,k).ne.HH(1,i,j,k).or.HH(2,i,j,k).ne.HH(2,i,j,k).or.HH(3,i,j,k).ne.HH(3,i,j,k).or.&
                   !   HH(3,i,j,k).ne.HH(3,i,j,k).or.HH(5,i,j,k).ne.HH(5,i,j,k)) then
                   !   print*, 'Flujos R H'
                   !endif
                else
                   call prim2f(fl, primF, 3)
                   call prim2f(fr, primB, 3)
                   HH(:,i,j,k) = (sr*fl(:)-sl*fr(:)+sl*sr*(UU(:,i,j,k+1)-UU(:,i,j,k)))/(sr-sl)
                   !if (HH(1,i,j,k).ne.HH(1,i,j,k).or.HH(2,i,j,k).ne.HH(2,i,j,k).or.HH(3,i,j,k).ne.HH(3,i,j,k).or.&
                   !   HH(3,i,j,k).ne.HH(3,i,j,k).or.HH(5,i,j,k).ne.HH(5,i,j,k)) then
                   !   print*, 'Flujos Mix H'
                   !endif
                endif
             enddo
          enddo
       enddo
       
    case(2)
       
       do i = nxmin+1, nx
          do j = nymin+1, ny
             do k = nzmin+1, nz
                
                primLL(:) = UPPrim(:,i-1,j,k)
                primL(:) = UPPrim(:,i,j,k)
                primR(:) = UPPrim(:,i+1,j,k)
                primRR(:) = UPPrim(:,i+2,j,k)
                call limiter(primLL,primL,primR,primRR,limtr)
                call wavespeed(primL, primR, sl, sr, 1)
                if (sl>=0) then
                   call prim2f(fff, primL, 1)
                   FF(:,i,j,k) = fff(:)
!                   if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
!                        .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                      !print*, 'Flujos L F'
!                      clF = clF + 1
!                   endif
                else if (sr<= 0) then
                   call prim2f(fff, primR, 1)
                   FF(:,i,j,k) = fff(:)
!                   if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
!                        .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                      !print*, 'Flujos R F'
!                      crF = crF + 1
!                   endif
                else
                   call prim2f(fl, primL, 1)
                   call prim2f(fr, primR, 1)
                   FF(:,i,j,k) = (sr*fl(:)-sl*fr(:)+sl*sr*(UU(:,i+1,j,k)-UU(:,i,j,k)))/(sr-sl)
!                   if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
!                        .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                      !print*, 'Flujos Mix F', primL(5), primR(5)
!                      cmF = cmF + 1
!                   endif
                endif
                
                primDD(:) = UPPrim(:,i,j-1,k)
                primD(:) = UPPrim(:,i,j,k) !! primD(:) = primL
                primU(:) = UPPrim(:,i,j+1,k)
                primUU(:) = UPPrim(:,i,j+2,k)
                call limiter(primDD,primD,primU,primUU,limtr)
                call wavespeed(primD, primU, sl, sr, 2)
                if (sl>=0) then
                   call prim2f(fff, primD, 2)
                   GG(:,i,j,k) = fff(:)
!                   if (GG(1,i,j,k).ne.GG(1,i,j,k).or.GG(2,i,j,k).ne.GG(2,i,j,k).or.GG(3,i,j,k).ne.GG(3,i,j,k).or.&
!                        GG(4,i,j,k).ne.GG(4,i,j,k).or.GG(5,i,j,k).ne.GG(5,i,j,k)) then
                      !print*, 'Flujos L G'
!                      clG = clG + 1
!                   endif
                else if (sr<=0) then
                   call prim2f(fff, primU, 2)
                   GG(:,i,j,k) = fff(:)
!                   if (GG(1,i,j,k).ne.GG(1,i,j,k).or.GG(2,i,j,k).ne.GG(2,i,j,k).or.GG(3,i,j,k).ne.GG(3,i,j,k).or.&
!                        GG(4,i,j,k).ne.GG(4,i,j,k).or.GG(5,i,j,k).ne.GG(5,i,j,k)) then
                      !print*, 'Flujos R G'
!                      crG = crG + 1
!                   endif
                else
                   call prim2f(fl, primD, 2)
                   call prim2f(fr, primU, 2)
                   GG(:,i,j,k) = (sr*fl(:)-sl*fr(:)+sl*sr*(UU(:,i,j+1,k)-UU(:,i,j,k)))/(sr-sl)
!                   if (GG(1,i,j,k).ne.GG(1,i,j,k).or.GG(2,i,j,k).ne.GG(2,i,j,k).or.GG(3,i,j,k).ne.GG(3,i,j,k).or.&
!                        GG(4,i,j,k).ne.GG(4,i,j,k).or.GG(5,i,j,k).ne.GG(5,i,j,k)) then
                      !print*, 'Flujos Mix G', primD(5), primU(5)
!                      cmG = cmG + 1
!                   endif
                endif
                
                primFF(:) = UPPrim(:,i,j,k-1)
                primF(:) = UPPrim(:,i,j,k) !! primF(:) = primL
                primB(:) = UPPrim(:,i,j,k+1)
                primBB(:) = UPPrim(:,i,j,k+2)
                call limiter(primFF,primF,primB,primBB,limtr)
                call wavespeed(primF, primB, sl, sr, 3)
                if (sl>=0) then
                   call prim2f(fff, primF, 3)
                   HH(:,i,j,k) = fff(:)
!                   if (HH(1,i,j,k).ne.HH(1,i,j,k).or.HH(2,i,j,k).ne.HH(2,i,j,k).or.HH(3,i,j,k).ne.HH(3,i,j,k).or.&
!                        HH(3,i,j,k).ne.HH(3,i,j,k).or.HH(5,i,j,k).ne.HH(5,i,j,k)) then
                      !print*, 'Flujos L H'
!                      clH = clH + 1
!                   endif
                else if (sr<=0) then
                   call prim2f(fff, primB, 3)
                   HH(:,i,j,k) = fff(:)
!                   if (HH(1,i,j,k).ne.HH(1,i,j,k).or.HH(2,i,j,k).ne.HH(2,i,j,k).or.HH(3,i,j,k).ne.HH(3,i,j,k).or.&
!                        HH(3,i,j,k).ne.HH(3,i,j,k).or.HH(5,i,j,k).ne.HH(5,i,j,k)) then
                      !print*, 'Flujos R H'
!                      crH = crH + 1
!                   endif
                else
                   call prim2f(fl, primF, 3)
                   call prim2f(fr, primB, 3)
                   HH(:,i,j,k) = (sr*fl(:)-sl*fr(:)+sl*sr*(UU(:,i,j,k+1)-UU(:,i,j,k)))/(sr-sl)
!                   if (HH(1,i,j,k).ne.HH(1,i,j,k).or.HH(2,i,j,k).ne.HH(2,i,j,k).or.HH(3,i,j,k).ne.HH(3,i,j,k).or.&
!                        HH(3,i,j,k).ne.HH(3,i,j,k).or.HH(5,i,j,k).ne.HH(5,i,j,k)) then
                      !print*, 'Flujos Mix H', primF(5), primB(5)
!                      cmH = cmH + 1
!                   endif
                endif
             enddo
          enddo
       enddo
!       print*, "clF = ", clF, "crF = ", crF, "cmF = ", cmF, "clG = ", clG, "crG = ", crG, "cmG = ", &
!            cmG, "clH = ", clH, "crH = ", crH, "cmH = ", cmH
    end select
    
  end subroutine HLLFluxes
  
  subroutine wavespeed(primL, primR, sl, sr, dim)
    
    use parameters, only: neq
    use HydroCore, only: csound
    
    implicit none
    
    integer, intent(in)                :: dim
    real, intent(in), dimension(neq)   :: primL, primR
    real, intent(out)                  :: sl, sr
    real                               :: csl, csr
    
    call csound(primL(5), primL(1), csl)
    call csound(primR(5), primR(1), csr)
    
    sl = min(primL(dim+1)-csl,primR(dim+1)-csr)
    sr = max(primL(dim+1)+csl,primR(dim+1)+csr)
    
  end subroutine wavespeed
  
  subroutine limiter (pll,pl,pr,prr,lim)
    
    use parameters, only: neq
    implicit none
    
    integer, intent(in) :: lim
    real, intent(in)    :: pll(neq)
    real, intent(in)    :: prr(neq)
    real, intent(inout) :: pl(neq)
    real, intent(inout) :: pr(neq)
    
    real                :: dl, dm, dr, al, ar
    integer             :: e
    
    do e = 1,neq
       dl = pl(e) - pll(e)
       dm = pr(e) - pl(e)
       dr = prr(e) - pr(e)
       al = average(dl, dm, lim)
       ar = average(dm, dr, lim)
       pl(e) = pl(e) + 0.5*al
       pr(e) = pr(e) - 0.5*ar
    end do
    
  contains
    
    real function average(a,b,opt)
      
      implicit none
      
      real, intent(in) :: a, b
      integer, intent(in) :: opt
      
      real :: s, c, d
      
      select case (opt)
         
      case (1) !LIMITER_NONE
         average = 0.5*(a+b)
         
      case (2) !LIMITER_VANLEER
         if (a*b.le.0.0) then
            average = 0.0
         else
            average = a*b*(a+b)/(a*a+b*b)
         end if
         
      case (3) !LIMITER_MINMOD
         s = sign(1.0,a)
         average = s*max(0.0, min(abs(a), s*b))
         
      case (4) !LIMITER_UMIST
         s = sign(1.0,a)
         c = 0.25*a + 0.75*b
         d = 0.75*a + 0.25*b
         average = min(2.0*abs(a), 2.0*s*b, s*c, s*d)
         average = s*max(0.0, average)
         
      case (5) !LIMITER_WOODWARD
         s = sign(1.0,a)
         c = 0.5*(a+b)
         average = min(2.0*abs(a), 2*s*b, s*c)
         average = s*max(0.0, average)
         
      case (6) !LIMITER_SUPERBEE
         s = sign(1.0,b)
         c = min(2.0*abs(b), s*a)
         d = min(abs(b),2.0*s*a)
         average = s*max(0.0,c,d)
         
      case default
         write(*,*) 'Nochingue compare'
         stop
         
      end select
      
    end function average
    
  end subroutine limiter
  
end module HLL
