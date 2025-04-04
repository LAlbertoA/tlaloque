  subroutine Godunov(UPPP,UU,FF,GG,HH,dt)

    use parameters, only: dx, dy, dz, nx, ny, nz, neq, nxmin, nymin, nzmin, nxmax, nymax, nzmax
#ifdef GRAV
    use sources
    use globals, only: PHI, UPrim, S, PHIP, PHIT
#endif
    implicit none

    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(out) :: UPPP
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(in) :: UU, FF, GG, HH
    real, intent(in)                :: dt
    real                            :: dtdx, dtdy, dtdz
    integer                         :: i,j,k,e!,png
#ifdef GRAV
    real, dimension(nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1) :: dens
    dens = UU(1,nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1)
    !dens = UU(1,nxmin+2:nxmax-2, nymin+2:nymax-2, nzmin+2:nzmax-2)

    !call poisson_solver1(dens, PHI, 3)
    call MultiGrid(dens,PHI)
    PHIT = PHIP + PHI
    call S_Solver(UPrim,PHIT,S)
#endif
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz
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

  subroutine HLLCFluxes(UU,FF,GG,HH,UPPrim,choice)

    use parameters, only: nxmin, nymin, nzmin, nx, ny, nz, neq, nxmax, nymax, nzmax, limtr
!    use HydroCore, only: prim2f
    
    implicit none
    
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(in) :: UU, UPPrim
    real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(out) :: FF, GG, HH
    integer, intent(in)               :: choice
    real, dimension(neq)              :: primLL,primL,primR,primRR,fff, fl, fr, primU, primB!, pl
    real, dimension(neq)              :: primDD, primUU, primBB, primFF, primD, primF, ustar!, pr
    real                              :: sl, sr, sst, rhostr
    integer                           :: i,j,k!,e,clF,crF,cmF,clG,crG,cmG,clH,crH,cmH
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
                call wavespeed(primL, primR, sl, sr, sst, 1)
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
                else if (sst >= 0) then
                   call prim2f(fl, primL, 1)
                   rhostr = primL(1)*(sl-primL(2))/(sl-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*sst
                   ustar(3) = rhostr*primL(3)
                   ustar(4) = rhostr*primL(4)
                   ustar(5) = rhostr*(UU(5,i,j,k)/primL(1) + (sst-primL(2))*(sst+primL(5)/(primL(1)*(sl-primL(2)))))
                   FF(:,i,j,k) = fl + sl*(ustar(:)-UU(:,i,j,k))
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos Mix F'
                   !endif
                else if (sst <= 0) then
                   call prim2f(fr, primR, 1)
                   rhostr = primR(1)*(sr-primR(2))/(sr-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*sst
                   ustar(3) = rhostr*primR(3)
                   ustar(4) = rhostr*primR(4)
                   ustar(5) = rhostr*(UU(5,i+1,j,k)/primR(1) + (sst-primR(2))*(sst+primR(5)/(primR(1)*(sr-primR(2)))))
                   FF(:,i,j,k) = fr + sr*(ustar(:)-UU(:,i+1,j,k))
                endif
                primD(:) = UPPrim(:,i,j,k)
                primU(:) = UPPrim(:,i,j+1,k)
                call wavespeed(primD, primU, sl, sr, sst, 2)
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
                else if (sst >= 0) then
                   call prim2f(fl, primD, 2)
                   rhostr = primD(1)*(sl-primD(3))/(sl-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primD(2)
                   ustar(3) = rhostr*sst
                   ustar(4) = rhostr*primD(4)
                   ustar(5) = rhostr*(UU(5,i,j,k)/primD(1) + (sst-primD(3))*(sst+primD(5)/(primD(1)*(sl-primD(3)))))
                   GG(:,i,j,k) = fl + sl*(ustar(:)-UU(:,i,j,k))
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos Mix F'
                   !endif
                else if (sst <= 0) then
                   call prim2f(fr, primU, 2)
                   rhostr = primU(1)*(sr-primU(3))/(sr-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primU(2)
                   ustar(3) = rhostr*sst
                   ustar(4) = rhostr*primU(4)
                   ustar(5) = rhostr*(UU(5,i,j+1,k)/primU(1) + (sst-primU(3))*(sst+primU(5)/(primU(1)*(sr-primU(3)))))
                   GG(:,i,j,k) = fr + sr*(ustar(:)-UU(:,i,j+1,k))
                endif
                primF(:) = UPPrim(:,i,j,k)
                primB(:) = UPPrim(:,i,j,k+1)
                call wavespeed(primF, primB, sl, sr, sst, 3)
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
                else if (sst >= 0) then
                   call prim2f(fl, primF, 3)
                   rhostr = primF(1)*(sl-primF(4))/(sl-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primF(2)
                   ustar(3) = rhostr*primF(3)
                   ustar(4) = rhostr*sst
                   ustar(5) = rhostr*(UU(5,i,j,k)/primF(1) + (sst-primF(4))*(sst+primF(5)/(primF(1)*(sl-primF(4)))))
                   HH(:,i,j,k) = fl + sl*(ustar(:)-UU(:,i,j,k))
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos Mix F'
                   !endif
                else if (sst <= 0) then
                   call prim2f(fr, primB, 3)
                   rhostr = primB(1)*(sr-primB(4))/(sr-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primB(2)
                   ustar(3) = rhostr*primB(3)
                   ustar(4) = rhostr*sst
                   ustar(5) = rhostr*(UU(5,i,j,k+1)/primB(1) + (sst-primB(4))*(sst+primB(5)/(primB(1)*(sr-primB(4)))))
                   HH(:,i,j,k) = fr + sr*(ustar(:)-UU(:,i,j,k+1))
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
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!! PPM TESTS !!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                if (i == 0) then
!                   call limiter(primLL,primL,primR,primRR,limtr)
!                if (i == nx) then
!                   call limiterppm(primLL,primL,primR,primRR,pr)
!                   call limiterppm(primL,primR,primRR,UPPrim(:,i-2,j,k),pl)
!                   do e = 1,neq
!                      if ((pl(e)-primL(e))*(primL(e)-pr(e)) <= 0) then
!                         pl(e) = primL(e)
!                         pr(e) = primL(e)
!                      endif
!                      if ((pl(e)-pr(e))*(primL(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pr(e) = 3.0*primL(e) - 2.0*pl(e)
!                      endif
!                      if (-(pl(e)-pr(e))*(primL(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pl(e) = 3.0*primL(e) - 2.0*pr(e)
!                      endif
!                   enddo
!                   primL(:) = pl(:)
!                   primR(:) = pr(:)
!                else
!                   call limiterppm(primLL,primL,primR,primRR,pr)
!                   call limiterppm(primL,primR,primRR,UPPrim(:,i+3,j,k),pl)
!                   do e = 1,neq
!                      if ((pl(e)-primL(e))*(primL(e)-pr(e)) <= 0) then
!                         pl(e) = primL(e)
!                         pr(e) = primL(e)
!                      endif
!                      if ((pl(e)-pr(e))*(primL(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pr(e) = 3.0*primL(e) - 2.0*pl(e)
!                      endif
!                      if (-(pl(e)-pr(e))*(primL(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pl(e) = 3.0*primL(e) - 2.0*pr(e)
!                      endif
!                   enddo
!                   primL(:) = pl(:)
!                   primR(:) = pr(:)
!                endif
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!! PPM TESTS !!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call limiter(primLL,primL,primR,primRR,limtr)
                call wavespeed(primL, primR, sl, sr, sst, 1)
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
                else if (sst >= 0) then
                   call prim2f(fl, primL, 1)
                   rhostr = primL(1)*(sl-primL(2))/(sl-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*sst
                   ustar(3) = rhostr*primL(3)
                   ustar(4) = rhostr*primL(4)
                   ustar(5) = rhostr*(UU(5,i,j,k)/primL(1) + (sst-primL(2))*(sst+primL(5)/(primL(1)*(sl-primL(2)))))
                   FF(:,i,j,k) = fl + sl*(ustar(:)-UU(:,i,j,k))
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
!                        .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                      !print*, 'Flujos Mix F', primL(5), primR(5)
!                      cmF = cmF + 1
!                   endif
                else if (sst <= 0) then
                   call prim2f(fr, primR, 1)
                   rhostr = primR(1)*(sr-primR(2))/(sr-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*sst
                   ustar(3) = rhostr*primR(3)
                   ustar(4) = rhostr*primR(4)
                   ustar(5) = rhostr*(UU(5,i+1,j,k)/primR(1) + (sst-primR(2))*(sst+primR(5)/(primR(1)*(sr-primR(2)))))
                   FF(:,i,j,k) = fr + sr*(ustar(:)-UU(:,i+1,j,k))
                endif
                primDD(:) = UPPrim(:,i,j-1,k)
                primD(:) = UPPrim(:,i,j,k) !! primD(:) = primL
                primU(:) = UPPrim(:,i,j+1,k)
                primUU(:) = UPPrim(:,i,j+2,k)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!! PPM TESTS !!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                if (j == 0) then
!                   call limiter(primDD,primD,primU,primUU,limtr)
!                if (j == ny) then
!                   call limiterppm(primDD,primD,primU,primUU,pr)
!                   call limiterppm(primD,primU,primUU,UPPrim(:,i,j-2,k),pl)
!                   do e = 1,neq
!                      if ((pl(e)-primD(e))*(primD(e)-pr(e)) <= 0) then
!                         pl(e) = primD(e)
!                         pr(e) = primD(e)
!                      endif
!                      if ((pl(e)-pr(e))*(primD(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pr(e) = 3.0*primD(e) - 2.0*pl(e)
!                      endif
!                      if (-(pl(e)-pr(e))*(primD(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pl(e) = 3.0*primD(e) - 2.0*pr(e)
!                      endif
!                   enddo
!                   primD(:) = pl(:)
!                   primU(:) = pr(:)
!                else
!                   call limiterppm(primDD,primD,primU,primUU,pr)
!                   call limiterppm(primD,primU,primUU,UPPrim(:,i,j+3,k),pl)
!                   do e = 1,neq
!                      if ((pl(e)-primD(e))*(primD(e)-pr(e)) <= 0) then
!                         pl(e) = primD(e)
!                         pr(e) = primD(e)
!                      endif
!                      if ((pl(e)-pr(e))*(primD(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pr(e) = 3.0*primD(e) - 2.0*pl(e)
!                      endif
!                      if (-(pl(e)-pr(e))*(primD(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pl(e) = 3.0*primD(e) - 2.0*pr(e)
!                      endif
!                   enddo
!                   primD(:) = pl(:)
!                   primU(:) = pr(:)
!                endif
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!! PPM TESTS !!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call limiter(primDD,primD,primU,primUU,limtr)
                call wavespeed(primD, primU, sl, sr, sst, 2)
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
                else if (sst >= 0) then
                   call prim2f(fl, primD, 2)
                   rhostr = primD(1)*(sl-primD(3))/(sl-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primD(2)
                   ustar(3) = rhostr*sst
                   ustar(4) = rhostr*primD(4)
                   ustar(5) = rhostr*(UU(5,i,j,k)/primD(1) + (sst-primD(3))*(sst+primD(5)/(primD(1)*(sl-primD(3)))))
                   GG(:,i,j,k) = fl + sl*(ustar(:)-UU(:,i,j,k))
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos Mix F'
                   !endif
                else if (sst <= 0) then
                   call prim2f(fr, primU, 2)
                   rhostr = primU(1)*(sr-primU(3))/(sr-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primU(2)
                   ustar(3) = rhostr*sst
                   ustar(4) = rhostr*primU(4)
                   ustar(5) = rhostr*(UU(5,i,j+1,k)/primU(1) + (sst-primU(3))*(sst+primU(5)/(primU(1)*(sr-primU(3)))))
                   GG(:,i,j,k) = fr + sr*(ustar(:)-UU(:,i,j+1,k))
                endif
                
                primFF(:) = UPPrim(:,i,j,k-1)
                primF(:) = UPPrim(:,i,j,k) !! primF(:) = primL
                primB(:) = UPPrim(:,i,j,k+1)
                primBB(:) = UPPrim(:,i,j,k+2)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!! PPM TESTS !!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                if (k == 0) then
!                   call limiter(primFF,primF,primB,primBB,limtr)
!                if (k == nz) then
!                   call limiterppm(primFF,primF,primB,primBB,pr)
!                   call limiterppm(primF,primB,primBB,UPPrim(:,i,j,k-2),pl)
!                   do e = 1,neq
!                      if ((pl(e)-primF(e))*(primF(e)-pr(e)) <= 0) then
!                         pl(e) = primF(e)
!                         pr(e) = primF(e)
!                      endif
!                      if ((pl(e)-pr(e))*(primF(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pr(e) = 3.0*primF(e) - 2.0*pl(e)
!                      endif
!                      if (-(pl(e)-pr(e))*(primF(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pl(e) = 3.0*primF(e) - 2.0*pr(e)
!                      endif
!                   enddo
!                   primF(:) = pl(:)
!                   primB(:) = pr(:)
!                else
!                   call limiterppm(primFF,primF,primB,primBB,pr)
!                   call limiterppm(primF,primB,primBB,UPPrim(:,i,j,k+3),pl)
!                   do e = 1,neq
!                      if ((pl(e)-primF(e))*(primF(e)-pr(e)) <= 0) then
!                         pl(e) = primF(e)
!                         pr(e) = primF(e)
!                      endif
!                      if ((pl(e)-pr(e))*(primF(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pr(e) = 3.0*primF(e) - 2.0*pl(e)
!                      endif
!                      if (-(pl(e)-pr(e))*(primF(e)-0.5*(pl(e)+pr(e))) > (pl(e)-pr(e))*(pl(e)-pr(e))/6.0) then
!                         pl(e) = 3.0*primF(e) - 2.0*pr(e)
!                      endif
!                   enddo
!                   primF(:) = pl(:)
!                   primB(:) = pr(:)
!                endif
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!! PPM TESTS !!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call limiter(primFF,primF,primB,primBB,limtr)
                call wavespeed(primF, primB, sl, sr, sst, 3)
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
                else if (sst >= 0) then
                   call prim2f(fl, primF, 3)
                   rhostr = primF(1)*(sl-primF(4))/(sl-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primF(2)
                   ustar(3) = rhostr*primF(3)
                   ustar(4) = rhostr*sst
                   ustar(5) = rhostr*(UU(5,i,j,k)/primF(1) + (sst-primF(4))*(sst+primF(5)/(primF(1)*(sl-primF(4)))))
                   HH(:,i,j,k) = fl + sl*(ustar(:)-UU(:,i,j,k))
                   !if (FF(1,i,j,k).ne.FF(1,i,j,k).or.FF(2,i,j,k).ne.FF(2,i,j,k).or.FF(3,i,j,k).ne.FF(3,i,j,k)&
                   !   .or.FF(4,i,j,k).ne.FF(4,i,j,k).or.FF(5,i,j,k).ne.FF(5,i,j,k)) then
                   !   print*, 'Flujos Mix F'
                   !endif
                else if (sst <= 0) then
                   call prim2f(fr, primB, 3)
                   rhostr = primB(1)*(sr-primB(4))/(sr-sst)
                   ustar(1) = rhostr
                   ustar(2) = rhostr*primB(2)
                   ustar(3) = rhostr*primB(3)
                   ustar(4) = rhostr*sst
                   ustar(5) = rhostr*(UU(5,i,j,k+1)/primB(1) + (sst-primB(4))*(sst+primB(5)/(primB(1)*(sr-primB(4)))))
                   HH(:,i,j,k) = fr + sr*(ustar(:)-UU(:,i,j,k+1))
                endif
             enddo
          enddo
       enddo
!       print*, "clF = ", clF, "crF = ", crF, "cmF = ", cmF, "clG = ", clG, "crG = ", crG, "cmG = ", &
!            cmG, "clH = ", clH, "crH = ", crH, "cmH = ", cmH
    end select
    
  end subroutine HLLCFluxes
  
  subroutine wavespeed(primL, primR, sl, sr, sst, dim)
    
    use parameters, only: neq
!    use HydroCore, only: csound
    
    implicit none
    
    integer, intent(in)                :: dim
    real, intent(in), dimension(neq)   :: primL, primR
    real, intent(out)                  :: sl, sr, sst
    real                               :: csl, csr
    
    call csound(primL(5), primL(1), csl)
    call csound(primR(5), primR(1), csr)
    
    sl = min(primL(dim+1)-csl,primR(dim+1)-csr)
    sr = max(primL(dim+1)+csl,primR(dim+1)+csr)
    sst = (primR(5)-primL(5) + primL(1)*primL(dim+1)*(sl-primL(dim+1)) - primR(1)*primR(dim+1)*(sr-primR(dim+1))) / &
         (primL(1)*(sl-primL(dim+1))-primR(1)*(sr-primR(dim+1)))

  end subroutine wavespeed

  subroutine limiterppm (pll,pl,pr,prr,p)
    
    use parameters, only: neq
    implicit none
    
    real, intent(in)    :: pll(neq)
    real, intent(in)    :: prr(neq)
    real, intent(in)    :: pl(neq)
    real, intent(in)    :: pr(neq)
    real, intent(out)   :: p(neq)
    
    real                :: dl, dm, dr, drr, dll, al, ar
    integer             :: e

    do e = 1,neq
       dl = pl(e) - pll(e)
       dm = pr(e) - pl(e)
       dr = prr(e) - pr(e)
       dll = pr(e) - pll(e)
       drr = prr(e) - pl(e)

       !Checamos por maximos o minimos
       
       if (dr*dm > 0) then
          ar = sign(1.0,drr)*min(0.5*abs(drr),2.0*abs(dm),2.0*abs(dr))
       else
          ar = 0.0
       endif

       if (dm*dl > 0) then
          al = sign(1.0,dll)*min(0.5*abs(dll),2.0*abs(dl),2.0*abs(dm))
       else
          al = 0.0
       endif

       p(e) = pl(e) + 0.5*dm - ((ar-al)/6.0)
    enddo

  end subroutine limiterppm
  
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
         average = min(2.0*abs(a), 2.0*s*b, s*c)
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
