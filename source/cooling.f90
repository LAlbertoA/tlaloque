module coolingh

  implicit none

contains

#ifdef COOL
  
  subroutine cooling()

    implicit none
    
    real :: maxloss
    
    ! Apply tabulated cooling to all local blocks
    
    call apply_cooling(maxloss)
    
  end subroutine cooling
  
  subroutine loaddata ()
    
    use parameters, only: cooling_file
    use globals, only: nptsT, cooltable, cool_Tmax, cool_Tmin
    implicit none
    
    integer :: i
    real :: a, b
    
    open (unit=99, file=cooling_file, status='old')
    
    read(99,*) nptsT
    allocate( cooltable(2,nptsT) )
    
    do i=1,nptsT
       read(99,*) a, b
       cooltable(1,i) = a
       cooltable(2,i) = -b
    end do
    close (unit=99)
    
    cool_Tmin = cooltable(1,1)
    cool_Tmax = cooltable(1,nptsT)
    
  end subroutine loaddata
  
  subroutine apply_cooling (maxloss)
    
    use parameters, only: nx, ny, nz, gamma, cooling_limit, mu0
    use globals, only: UPrim, U, cool_Tmin, dt
    use constants, only: AMU
    implicit none
    
    real, intent(out) :: maxloss
    
    real, parameter :: T_floor = 10.0
    
    real :: new_temp, temp, vel2, ETH, EK, aloss
    real :: frac_loss, numdens, ce, cool_factor
    integer :: i, j, k
    
    maxloss = 0.0
    
    do i=1,nx
       do j=1,ny
          do k=1,nz
             
             ! Calculate temperature of this cell
             call calcTemp (UPrim(:,i,j,k), temp)
             ! Cooling not applied below cool_Tmin
             if (temp.ge.cool_Tmin) then
                
                ! Code units
                ETH = UPrim(5,i,j,k)/(gamma-1)
                vel2 = UPrim(2,i,j,k)**2 + UPrim(3,i,j,k)**2 + UPrim(4,i,j,k)**2
                EK = 0.5 * UPrim(1,i,j,k) * vel2
                
                ! Interpolate cooling coefficient from table
                call find_coef (temp, aloss)
                
                ! Calculate radiative loss and new temperature
                numdens = UPrim(1,i,j,k)/(mu0*AMU)  ! cgs, gas ionized
                ce = aloss*(numdens**2)/ETH  ! cgs
                cool_factor = exp(-ce*(dt))
                frac_loss = 1.0-cool_factor
                
                ! DEBUG
                !          write(logu,*) localBlocks(bIndx), i, j, k
                !          write(logu,*) log10(temp), aloss
                !          write(logu,*) numdens, frac_loss
                
                ! Record maximum cooling for this block before limiting
                maxloss = max(maxloss, frac_loss)
                
                ! Limit cool_factor directly, if needed
                if (cool_factor.lt.1.0-cooling_limit) then
                   cool_factor = 1.0-cooling_limit
                end if
                
                ! Impose a temperature floor by adjusting cool_factor, if needed
                ! Set to 10 K by default
                new_temp = temp * cool_factor
                if (new_temp.lt.T_floor) then
                   new_temp = T_floor
                   cool_factor = T_floor / temp
                end if
                
                ! Update pressure and total energy
                UPrim(5,i,j,k) = UPrim(5,i,j,k) * cool_factor
                ETH = UPrim(5,i,j,k)/(gamma-1)
                U(5,i,j,k) = EK + ETH
                
             end if
             
          end do
       end do
    end do
    
  end subroutine apply_cooling
  
  subroutine find_coef (temp, coef)
    
    use globals, only: cooltable, nptsT
    implicit none
    
    real, intent(in) :: temp
    real, intent(out) :: coef
    
    integer :: i
    real :: T0, T1, C0, C1
    
    ! For T in range of the table, we do linear interpolation
    do i=2,nptsT
       if (cooltable(1,i).gt.temp) then
          T0 = log10(cooltable(1,i-1))
          C0 = cooltable(2,i-1)
          T1 = log10(cooltable(1,i))
          C1 = cooltable(2,i)
          coef = 10**( (C1-C0)/(T1-T0)*(log10(temp)-T0) + C0 )
          return
       end if
    end do
    
  end subroutine find_coef
  
  subroutine calcTemp (prims, temp)
    
    use parameters, only: mu0, neq
    use constants, only: KB, AMU
    implicit none
    
    real, intent(in)  :: prims(neq)
    real, intent(out) :: temp
    
    temp = prims(5)*mu0*AMU/(prims(1)*KB)
    
    return
    
  end subroutine calcTemp

#endif
  
end module coolingh
