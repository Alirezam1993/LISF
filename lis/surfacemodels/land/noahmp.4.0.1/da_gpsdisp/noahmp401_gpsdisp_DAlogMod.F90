!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! Oct 2024: Alireza Moghaddasi; Created for GPS displacement DA with daily logging
!

module noahmp401_gpsdisp_DAlogMod
  
  use LIS_constantsMod,  only : LIS_CONST_RHOFW
  use ESMF
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------
  public :: noahmp401_gpsdisp_DAlog
  public :: initialize_gps_structures
!-----------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------
  public :: NOAHMPgpspred_struc
!EOP

  type, public :: NOAHMPgpspred_dec
     real, allocatable :: daily_tws(:)  ! Store daily TWS for current day
  end type NOAHMPgpspred_dec
  
  type (NOAHMPgpspred_dec), allocatable :: NOAHMPgpspred_struc(:)
  
  logical :: structures_initialized = .false.
  
contains 
  
  subroutine noahmp401_gpsdisp_DAlog(n)
    
    ! USES:
    use LIS_coreMod, only : LIS_rc, LIS_surface
    use LIS_timeMgrMod
    use noahmp401_lsmMod
    use LIS_logMod, only : LIS_logunit, LIS_verify
      
    ! ARGUMENTS:  
    integer, intent(in) :: n 
      
    ! DESCRIPTION:
    ! Calculates total water storage DAILY for GPS displacement data assimilation.
    ! Simplified approach following GRACE DA pattern.
    
    integer :: t          
    real :: current_tws

    ! Only log during forward mode (not during increment application)
    if (LIS_rc%DAincrMode(n).eq.0) then
       
       ! Initialize structures if needed
       if (.not. structures_initialized) then
          call initialize_gps_structures(n)
       endif
       
       ! Log TWS at specific times (following GRACE pattern but more frequent)
       ! Log daily at 12:00:00 to match GRACE timing structure
       if((LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then
          write(LIS_logunit,*) '[INFO] GPS-DA: Daily TWS logging for day:', LIS_rc%da, 'hour:', LIS_rc%hr

          ! Calculate and store TWS for each tile
          do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
             current_tws = &
                  NOAHMP401_struc(n)%noahmp401(t)%sneqv   +         &
                  (NOAHMP401_struc(n)%noahmp401(t)%canliq  +         &
                  NOAHMP401_struc(n)%noahmp401(t)%canice) +         &
                  (NOAHMP401_struc(n)%noahmp401(t)%smc(1)  *         &
                  NOAHMP401_struc(n)%sldpth(1)*LIS_CONST_RHOFW)        +         &
                  (NOAHMP401_struc(n)%noahmp401(t)%smc(2)  *         &
                  NOAHMP401_struc(n)%sldpth(2)*LIS_CONST_RHOFW)        +         &
                  (NOAHMP401_struc(n)%noahmp401(t)%smc(3)  *         &
                  NOAHMP401_struc(n)%sldpth(3)*LIS_CONST_RHOFW)        +         &
                  (NOAHMP401_struc(n)%noahmp401(t)%smc(4)  *         &
                  NOAHMP401_struc(n)%sldpth(4)*LIS_CONST_RHOFW)        +         &
                  NOAHMP401_struc(n)%noahmp401(t)%wa
             
             ! Store daily TWS directly (no averaging needed for GPS DA)
             NOAHMPgpspred_struc(n)%daily_tws(t) = current_tws
          enddo
          
          write(LIS_logunit,*) '[INFO] GPS-DA: Daily TWS logging completed'
       endif
    endif

  end subroutine noahmp401_gpsdisp_DAlog
  
  subroutine initialize_gps_structures(n)
    ! USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_logMod, only : LIS_logunit
    
    ! ARGUMENTS:
    integer, intent(in) :: n
    
    ! DESCRIPTION:
    ! Initialize GPS DA structures (simplified)
    
    integer :: m
    
    if (.not. structures_initialized) then
      ! Initialize prediction structures for ALL nests
      if (.not.allocated(NOAHMPgpspred_struc)) then 
         allocate(NOAHMPgpspred_struc(LIS_rc%nnest))
         do m = 1, LIS_rc%nnest
            allocate(NOAHMPgpspred_struc(m)%daily_tws(&
                 LIS_rc%npatch(m,LIS_rc%lsm_index)))
            NOAHMPgpspred_struc(m)%daily_tws = 0.0
         end do
         write(LIS_logunit,*) '[INFO] GPS-DA: Allocated prediction structures'
      endif
      
      structures_initialized = .true.
    endif
    
  end subroutine initialize_gps_structures
  
end module noahmp401_gpsdisp_DAlogMod 