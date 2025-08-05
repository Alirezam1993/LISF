!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: GPSdispObs_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle GPS vertical displacement observations
!   
! !REVISION HISTORY: 
!  27Feb05    Sujay Kumar;   Initial Specification
!  14Aug20    Gaohong Yin;   Changed for GPS displacement DA
!  24Oct24    Alireza Moghaddasi; Revised to align with GRACE code
! 
module GPSdispObs_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_perturbMod
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GPSdispObs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GPS_struc
  PUBLIC :: GPS_observation

  ! Define the GPS_observation type
  type, public :: GPS_observation
     real    :: value = -9999.0  ! Observed displacement value
     real    :: error = -9999.0  ! Error/uncertainty in observation
     logical :: assim = .false.  ! Whether to assimilate this observation
  end type GPS_observation

  type, public :: GPS_dec
     integer           :: mo
     integer           :: alarmhr
     integer           :: useDistErr
  end type GPS_dec

  type(GPS_dec), allocatable :: GPS_struc(:)
  
  ! Observations array accessible to other modules
  type(GPS_observation), allocatable, public :: Observations(:)

contains

!BOP
! 
! !ROUTINE: GPSdispObs_setup
! \label{GPSdispObs_setup}
! 
! !INTERFACE: 
  subroutine GPSdispObs_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_DAobservationsMod
    use LIS_perturbMod

    implicit none 

! !ARGUMENTS: 
    integer, intent(IN)    ::  k  ! Dataset index
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION:
!   
!   This routine completes the runtime initializations and 
!   creation of data structures required for GPS displacement assimilation.
!  
!   The arguments are: 
!   \begin{description}
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n
    integer                ::  ftn
    integer                ::  i
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) :: GPSdispobsdir
    character*100          ::  temp

    real, allocatable      :: ssdev(:)
    character*1            :: vid(2)
    character*40, allocatable :: vname(:)
    real, allocatable      :: varmin(:)
    real, allocatable      :: varmax(:)
    type(pert_dec_type)    :: obs_pert
    real, pointer          :: obs_temp(:,:)

    write(LIS_logunit,*) '[GPS_DA] Setting up GPS data assimilation'

        allocate(GPS_struc(LIS_rc%nnest))

    ! Set up array specifications
    call ESMF_ArraySpecSet(intarrspec, rank=1, typekind=ESMF_TYPEKIND_I4, rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: GPSdispObs_setup')

    call ESMF_ArraySpecSet(realarrspec, rank=1, typekind=ESMF_TYPEKIND_R4, rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: GPSdispObs_setup')

    call ESMF_ArraySpecSet(pertArrSpec, rank=2, typekind=ESMF_TYPEKIND_R4, rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: GPSdispObs_setup')

    ! Get GPS data directory
    call ESMF_ConfigFindLabel(LIS_config, "GPS Displacement data directory:", rc=status)
    call LIS_verify(status, 'GPS Displacement data directory label not found')

    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, GPSdispobsdir, rc=status)
        call LIS_verify(status, 'GPS Displacement data directory not defined')
        call ESMF_AttributeSet(OBS_State(n), "Data Directory", GPSdispobsdir, rc=status)
        call LIS_verify(status, 'Error setting Data Directory attribute')
    enddo

    ! Get configuration for GPS measurement error
    call ESMF_ConfigFindLabel(LIS_config, "GPS use reported measurement error values:", rc=status)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, GPS_struc(n)%useDistErr, rc=status)
        call LIS_verify(status, 'GPS use reported measurement error values: not defined')
    enddo

    ! Set initial attributes in OBS_State
    do n=1, LIS_rc%nnest
        call ESMF_AttributeSet(OBS_State(n), "Data Update Status", .false., rc=status)
        call LIS_verify(status, 'Error setting Data Update Status')

        call ESMF_AttributeSet(OBS_State(n), "Data Update Time", -99.0, rc=status)
        call LIS_verify(status, 'Error setting Data Update Time')

        call ESMF_AttributeSet(OBS_State(n), "Data Assimilate Status", .false., rc=status)
        call LIS_verify(status, 'Error setting Data Assimilate Status')

        call ESMF_AttributeSet(OBS_State(n), "Number Of Observations", LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error setting Number Of Observations')
    enddo

    write(LIS_logunit,*) '[INFO] Read GPS displacement data specifications'

    !----------------------------------------------------------------------------
    !   Create the array containers that will contain the observations and
    !   the perturbations. For the GPS displacement case, it is assumed that the 
    !   observations are in the grid space, similar to GRACE. 
    !----------------------------------------------------------------------------

    do n=1, LIS_rc%nnest
        write(unit=temp, fmt='(i2.2)') 1
        read(unit=temp, fmt='(2a1)') vid

        obsField(n) = ESMF_FieldCreate( &
             arrayspec=realarrspec, grid=LIS_obsvecGrid(n,k), &
             name="Observation"//vid(1)//vid(2), rc=status)
        call LIS_verify(status, 'Error ESMF_FieldCreate: GPSdispObs_setup')

        ! Perturbations State
        write(LIS_logunit,*) '[INFO] Opening attributes for observations ', &
             trim(LIS_rc%obsattribfile(k))
        ftn = LIS_getNextUnitNumber()
        open(ftn, file=trim(LIS_rc%obsattribfile(k)), status='old')
        read(ftn,*)
        read(ftn,*) LIS_rc%nobtypes(k)
        read(ftn,*)
    
        allocate(vname(LIS_rc%nobtypes(k)))
        allocate(varmax(LIS_rc%nobtypes(k)))
        allocate(varmin(LIS_rc%nobtypes(k)))

        do i=1, LIS_rc%nobtypes(k)
            read(ftn, fmt='(a40)') vname(i)
            read(ftn,*) varmin(i), varmax(i)
            write(LIS_logunit,*) '[INFO] ', vname(i), varmin(i), varmax(i)
        enddo
        call LIS_releaseUnitNumber(ftn)

        allocate(ssdev(LIS_rc%obs_ngrid(k)))

        if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
            allocate(obs_pert%vname(1))
            allocate(obs_pert%perttype(1))
            allocate(obs_pert%ssdev(1))
            allocate(obs_pert%stdmax(1))
            allocate(obs_pert%zeromean(1))
            allocate(obs_pert%tcorr(1))
            allocate(obs_pert%xcorr(1))
            allocate(obs_pert%ycorr(1))
            allocate(obs_pert%ccorr(1,1))

            call LIS_readPertAttributes(1, LIS_rc%obspertAttribfile(k), obs_pert)

            ssdev = obs_pert%ssdev(1)

            pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec, &
                 grid=LIS_obsEnsOnGrid(n,k), &
                 name="Observation"//vid(1)//vid(2), &
                 rc=status)
            call LIS_verify(status, 'Error ESMF_FieldCreate: GPSdispObs_setup')

            ! Initializing the perturbations to be zero 
            call ESMF_FieldGet(pertField(n), localDE=0, farrayPtr=obs_temp, rc=status)
            call LIS_verify(status, 'Error ESMF_FieldGet: GPSdispObs_setup')
            obs_temp(:,:) = 0

            call ESMF_AttributeSet(pertField(n), "Perturbation Type", &
                 obs_pert%perttype(1), rc=status)
            call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')

            if(LIS_rc%obs_ngrid(k).gt.0) then 
                call ESMF_AttributeSet(pertField(n), "Standard Deviation", &
                     ssdev, itemCount=LIS_rc%obs_ngrid(k), rc=status)
                call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')
            endif

            call ESMF_AttributeSet(pertField(n), "Std Normal Max", &
                 obs_pert%stdmax(1), rc=status)
            call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')

            call ESMF_AttributeSet(pertField(n), "Ensure Zero Mean", &
                 obs_pert%zeromean(1), rc=status)
            call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')

            call ESMF_AttributeSet(pertField(n), "Temporal Correlation Scale", &
                 obs_pert%tcorr(1), rc=status)
            call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')

            call ESMF_AttributeSet(pertField(n), "X Correlation Scale", &
                 obs_pert%xcorr(1), rc=status)
            call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')

            call ESMF_AttributeSet(pertField(n), "Y Correlation Scale", &
                 obs_pert%ycorr(1), rc=status)
            call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')

            call ESMF_AttributeSet(pertField(n), "Cross Correlation Strength", &
                 obs_pert%ccorr(1,:), itemCount=1, rc=status)
            call LIS_verify(status, 'Error ESMF_AttributeSet: GPSdispObs_setup')
        endif

        deallocate(vname)
        deallocate(varmax)
        deallocate(varmin)
        deallocate(ssdev)
    enddo

    write(LIS_logunit,*) '[INFO] Created the States to hold the GPS displacement observations data'

!--------------------------------------------------------------------------------
! The data will be read and kept in memory at 3z, similar to GRACE data
!--------------------------------------------------------------------------------    
    do n=1, LIS_rc%nnest
        GPS_struc(n)%alarmHr = 12  ! GPS DA at noon to match TWS logging

        call ESMF_StateAdd(OBS_State(n), (/obsField(n)/), rc=status)
        call LIS_verify(status, 'Error ESMF_StateAdd: GPSdispObs_setup')

        call ESMF_StateAdd(OBS_Pert_State(n), (/pertField(n)/), rc=status)
        call LIS_verify(status, 'Error ESMF_StateAdd: GPSdispObs_setup')  
    enddo
    
    ! Initialize Observations array
    if (LIS_rc%obs_ngrid(k) > 0) then
        if (allocated(Observations)) deallocate(Observations)
        allocate(Observations(LIS_rc%obs_ngrid(k)))
        Observations(:)%value = -9999.0
        Observations(:)%error = -9999.0
        Observations(:)%assim = .false.
    endif
    
    write(LIS_logunit,*) '[INFO] GPS data assimilation setup complete'

  end subroutine GPSdispObs_setup

end module GPSdispObs_Mod
