program HydroCode

    use parameters
    use globals
    use constants
#ifdef COOL
    use coolingh
#endif
    use winds
    use pgravity
    
    implicit none 

#ifdef GRAV
    real, dimension(nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1)  :: dens
#endif
    !!! ------------------ main ----------------------------------------

    it = 0
    nout = 0
    itprint = 0
    tout = 0.0
    dt = 0.0
    t = 0.0
    dth = 0.0

    call initmain(nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq)
    
    call initflow()

    if (rank == 0) then
       print*, 'Condiciones iniciales listas'
    endif
    
    call borders(U,1)

    if (rank == 0) then
       print*, 'Fronteras listas'
    endif

    call calc_prims(UPrim,U)

    if (rank == 0) then
       print*, 'Primitivas listas'
    endif
    
#ifdef GRAV    
    dens = U(1,nxmin+1:nxmax-1, nymin+1:nymax-1, nzmin+1:nzmax-1)
    call Multigrid(dens,PHI)
    if (rank == 0) then
       print*, 'Potencial listo'
    endif
#endif    
    if (rank == 0) then
        print*, 'Inicia ciclo principal'
    endif

    do while (t<tfin)

       if (rank == 0) then
          print *, 'dt = ', dt, ' t = ', t, ' it = ', it
          if (logged.eqv..true.) then
             write(logu,*) 'dt = ', dt, ' t = ', t, ' it = ', it
          endif
       endif
       
       if (t>=tout) then
          call output(nout, tout)
       endif
              
       call time_step(UPrim, dt)
       !        print*, '1'
       call HLLCFluxes(U,F,G,H,UPrim,1)
       !        print*, '2'
       dth = dt/2.0
       !
       call Godunov(UPP,U,F,G,H,dth)
       !        print*, '3'
       call borders(UPP,2)
       !        print*, '4'
       call calc_prims(UPrim,UPP)
       !        print*, '5'
       call HLLCFluxes(UPP,F,G,H,UPrim,2)
       !        print*, '6'
       call Godunov(UP,U,F,G,H,dt)
       !        print*, '7'
       call borders(UP,1)
       !        print*, '8'
       call step(t,it,dt)
       !        print*, '9'
       call calc_prims(UPrim,U)
       !        print*, '10'
#ifdef COOL
       call cooling()
       !        print*, '11'
#endif
    enddo
    !
    call output(nout, tout)
    !
    if (rank==0) then
       print*, "Done"
       print*, 'dt = ', dt, ' t = ', t, ' it = ', it
       if (logged.eqv..true.) then
          close(logu)
       endif
    endif
#ifdef MPIP
    call mpi_finalize(err)
#endif

    !!! ------------------ main ----------------------------------------

end program
