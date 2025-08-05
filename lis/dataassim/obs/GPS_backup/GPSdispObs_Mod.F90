!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
!  29Sep17    Yonghwan Kwon;  Modified for DTB observations
!  14Aug20    Gaohong Yin; Changed for GPS displacement DA 
!  24Oct24    Alireza Moghaddasi; Revised to align with GRACE code
! 
module GPSdispObs_Mod
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  ! Remove or comment this import to avoid conflict with our own GPS_observation type
  ! use enksgps_types                ! Contains obs_type definition
  use LIS_perturbMod
!EOP
  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GPSdispObs_setup
  ! Add these new functions to public interface
  PUBLIC :: is_observations_valid
  PUBLIC :: safely_allocate_observations
  PUBLIC :: safely_deallocate_observations
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GPS_struc
  ! Define our Observations type as PUBLIC
  PUBLIC :: GPS_observation
  ! Make the Observations array variable public
  PUBLIC :: Observations

  type, public :: GPS_dec
     integer           :: mo
     integer           :: alarmhr
     integer           :: useDistErr
  end type GPS_dec

  type(GPS_dec), allocatable :: GPS_struc(:)

  ! Define a structure to store GPS displacement observations
  type :: GPS_observation
     real :: value       ! The observed displacement value
     real :: error       ! Error/uncertainty associated with this observation
     logical :: assim    ! Flag indicating if this observation should be assimilated
  end type GPS_observation

  ! Define the public observations array for use in data assimilation
  type(GPS_observation), allocatable :: Observations(:)

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
    integer                ::  Nobjs, N_obs_size, Nobs, NobsCheck, gridNum, coordCount
    character(len=20)      :: nest_str, ngrid_str  ! String conversion buffers
    character(len=10)      :: GRIDTYPE
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
    logical                :: file_exists
    logical                :: alloc_success ! For safe allocation

    ! Initial debug statements
    write(LIS_logunit,*) '[GPS_DA] ENTERING SETUP ROUTINE'
    write(LIS_logunit,*) '[GPS_DA] Total nests configured (nnest):', LIS_rc%nnest
    write(LIS_logunit,*) '[GPS_DA] obs_ngrid array size:', size(LIS_rc%obs_ngrid)
    write(LIS_logunit,*) '[GPS_DA] obs_ngrid values:', LIS_rc%obs_ngrid
    write(LIS_logunit,*) '[GPS_DA] Current dataset index (k):', k
    if (k > size(LIS_rc%obs_ngrid)) then
        write(LIS_logunit,*) '[ERR] Dataset index k exceeds obs_ngrid size!'
        call LIS_endrun()
    endif

    if (size(LIS_rc%obs_ngrid) < LIS_rc%nnest) then
        write(LIS_logunit,*) '[ERR] obs_ngrid size mismatch! nnest:',LIS_rc%nnest,&
                             'obs_ngrid size:',size(LIS_rc%obs_ngrid)
        call LIS_endrun()
    endif

    ! Validate nest configuration
    if (LIS_rc%nnest <= 0) then
        write(LIS_logunit,*) '[ERR] Invalid number of nests:', LIS_rc%nnest
        call LIS_endrun()
    endif

    if (size(LIS_rc%obs_ngrid) /= LIS_rc%nnest) then
        write(LIS_logunit,*) '[ERR] Config mismatch: nnest=', LIS_rc%nnest, &
                             ' vs obs_ngrid size=', size(LIS_rc%obs_ngrid)
        call LIS_endrun()
    endif

    do n=1, LIS_rc%nnest
        if (LIS_rc%obs_ngrid(n) <= 0) then
            write(LIS_logunit,*) '[ERR] Invalid obs_ngrid(',n,')=', LIS_rc%obs_ngrid(n)
            call LIS_endrun()
        endif
    enddo

    ! Attempt to retrieve GRIDTYPE attribute
    do n=1, LIS_rc%nnest
        call ESMF_AttributeGet(OBS_State(n), "GRIDTYPE", GRIDTYPE, rc=status)
        if (status /= ESMF_SUCCESS) then
            GRIDTYPE = '1D'
            write(nest_str,'(I0)') n
            call ESMF_AttributeSet(OBS_State(n), "GRIDTYPE", GRIDTYPE, rc=status)
            call LIS_verify(status,'Failed setting GRIDTYPE for nest='//trim(nest_str))
        endif
    enddo

    ! Add before validation check
    write(LIS_logunit,*) '[GPS_DA] Configuration Parameters:'
    write(LIS_logunit,*) '         obs_ngrid =', LIS_rc%obs_ngrid(k)
    write(LIS_logunit,*) '         gnc        =', LIS_rc%gnc(k)
    write(LIS_logunit,*) '         gnr        =', LIS_rc%gnr(k)

    ! Add after setting obs_ngrid
    write(LIS_logunit,*) '[GPS_DA] Observation points:', LIS_rc%obs_ngrid(k)

    ! Allocate GPS structure
    if (.not. allocated(GPS_struc)) then
        allocate(GPS_struc(LIS_rc%nnest))
        write(LIS_logunit,*) '[GPS_DA] Allocated GPS_struc for', LIS_rc%nnest, 'nests'
        if (LIS_rc%nnest <= 0) then
            write(LIS_logunit,*) '[ERR] Invalid number of nests:', LIS_rc%nnest
            call LIS_endrun()
        endif
    endif

    ! Set up array specifications
    call ESMF_ArraySpecSet(intarrspec, rank=1, typekind=ESMF_TYPEKIND_I4, rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: GPSdispObs_setup')

    call ESMF_ArraySpecSet(realarrspec, rank=1, typekind=ESMF_TYPEKIND_R4, rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: GPSdispObs_setup')

    call ESMF_ArraySpecSet(pertArrSpec, rank=2, typekind=ESMF_TYPEKIND_R4, rc=status)
    call LIS_verify(status, 'Error ESMF_ArraySpecSet: GPSdispObs_setup')

    ! Don't try to initialize pertField here - we'll handle it when it's used
    ! Just note that it's uninitialized at this point
    write(LIS_logunit,*) '[DEBUG] Note: pertField is intentionally uninitialized until needed'

    ! Get GPS data directory
    call ESMF_ConfigFindLabel(LIS_config, "GPS Displacement data directory:", rc=status)
    call LIS_verify(status, 'GPS Displacement data directory label not found')

    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, GPSdispobsdir, rc=status)
        call LIS_verify(status, 'GPS Displacement data directory not defined')
        call ESMF_AttributeSet(OBS_State(n), "Data Directory", GPSdispobsdir, rc=status)
        call LIS_verify(status, 'Error setting Data Directory attribute')
        write(LIS_logunit,*) '[DEBUG] Set Data Directory to:', trim(GPSdispobsdir)
    enddo

    ! Set initial attributes in OBS_State
    do n=1, LIS_rc%nnest
        call ESMF_AttributeSet(OBS_State(n), "Data Update Status", .false., rc=status)
        call LIS_verify(status, 'Error setting Data Update Status')
        write(LIS_logunit,*) '[DEBUG] Set Data Update Status to: .false.'

        call ESMF_AttributeSet(OBS_State(n), "Data Update Time", -99.0, rc=status)
        call LIS_verify(status, 'Error setting Data Update Time')
        write(LIS_logunit,*) '[DEBUG] Set Data Update Time to: -99.0'

        call ESMF_AttributeSet(OBS_State(n), "Data Assimilate Status", .false., rc=status)
        call LIS_verify(status, 'Error setting Data Assimilate Status')
        write(LIS_logunit,*) '[DEBUG] Set Data Assimilate Status to: .false.'

        ! Calculate Number of Observations
        Nobjs = 1  ! Example value; ensure this aligns with your logic
        N_obs_size = 46  ! Example value; ensure this aligns with your logic
        Nobs = Nobjs * N_obs_size

        ! 1) Set the actual number of observations
        write(nest_str,'(I0)') n
        write(ngrid_str,'(I0)') LIS_rc%obs_ngrid(n)
        call ESMF_AttributeSet(OBS_State(n), "GPS_NumObs", LIS_rc%obs_ngrid(n), rc=status)
        call LIS_verify(status,'GPS_NumObs failed nest='//trim(nest_str)//&
                             ' value='//trim(ngrid_str))

        ! Initialize Observations array properly using our safe function
        call safely_allocate_observations(LIS_rc%obs_ngrid(n))

        ! Set GPS structure
        GPS_struc(n)%alarmhr = 3
        write(LIS_logunit,*) '[DEBUG] Set GPS alarm hour to:', GPS_struc(n)%alarmhr
    enddo

    write(LIS_logunit,*) '[INFO] Created the States to hold the GPS observations data'

    ! Check to make sure OBS_State fields are initialized
    write(LIS_logunit,*) '[DEBUG] Checking OBS_State initialization'
    do n=1, LIS_rc%nnest
        write(nest_str,'(I0)') n
        ! Ensure OBS_State has required attributes
        call ESMF_AttributeSet(OBS_State(n), "GRIDTYPE", "1D", rc=status)
        call LIS_verify(status, 'Failed setting GRIDTYPE attribute in OBS_State')
    enddo
    
    ! This code is premature - pertField is not created yet, so commenting it out
    ! Ensure attributes are set on PERTURBATION fields
    !do n=1, LIS_rc%nnest
    !    write(nest_str,'(I0)') n
    !    if (trim(LIS_rc%perturb_obs(n)) /= "none") then
    !        call ESMF_AttributeSet(pertField(n),"Grid Number",LIS_rc%obs_ngrid(n),itemCount=LIS_rc%obs_ngrid(n),rc=status)
    !        call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
    !        write(LIS_logunit,*) '[DEBUG] Set Grid Number:', LIS_rc%obs_ngrid(n)
    !    endif
    !enddo

    ! Extra debugging for safer execution
    write(LIS_logunit,*) '[DEBUG] About to write grid configuration'

    ! Print basic grid info to help diagnose issues
    write(LIS_logunit,*) '[GPS_DA] Grid configuration:'
    write(LIS_logunit,*) '         Columns:', LIS_rc%gnc(k)
    write(LIS_logunit,*) '         Rows:   ', LIS_rc%gnr(k)
    write(LIS_logunit,*) '         Total points:', LIS_rc%obs_ngrid(k)

    ! Add more debug info and safety check
    write(LIS_logunit,*) '[DEBUG] About to check Observations array'

    ! Fix the observation processing loop with added safety checks
    if (allocated(Observations)) then
        write(LIS_logunit,*) '[DEBUG] Observations array is allocated with size:', size(Observations)
        
        ! Additional safety check for array bounds
        if (LIS_rc%obs_ngrid(k) > 0 .and. LIS_rc%obs_ngrid(k) <= size(Observations)) then
            do i = 1, LIS_rc%obs_ngrid(k)
                if (Observations(i)%assim) then
                    write(LIS_logunit,*) '[GPS_DA] Processing observation at point:', i
                endif
            enddo
        else
            write(LIS_logunit,*) '[WARN] Cannot process observations: invalid size comparison', &
                                 LIS_rc%obs_ngrid(k), 'vs', size(Observations)
        endif
    else
        write(LIS_logunit,*) '[WARN] Observations array not allocated'
    endif

    ! Add debug context to critical operations
    do n=1, LIS_rc%nnest
        write(LIS_logunit,*) '[DEBUG] Nest:',n,' Configuration:'
        write(LIS_logunit,*) '  obs_ngrid:', LIS_rc%obs_ngrid(n)
        write(LIS_logunit,*) '  gnc:', LIS_rc%gnc(n)
        write(LIS_logunit,*) '  gnr:', LIS_rc%gnr(n)
        if (allocated(LIS_rc%perturb_obs)) then
            write(LIS_logunit,*) '  perturb_obs:', trim(LIS_rc%perturb_obs(n))
        else
            write(LIS_logunit,*) '  perturb_obs: NOT ALLOCATED'
        endif
    enddo

    write(LIS_logunit,*) '[INFO] Read GPS displacement data specifications'

    !----------------------------------------------------------------------------
    !   Create the array containers that will contain the observations and
    !   the perturbations.
    !----------------------------------------------------------------------------
    do n=1, LIS_rc%nnest
        write(unit=temp, fmt='(i2.2)') 1
        read(unit=temp, fmt='(2a1)') vid

        ! Debug output
        write(LIS_logunit,*) '[DEBUG] Creating field for nest:', n
        write(LIS_logunit,*) '[DEBUG] Grid size:', LIS_rc%obs_ngrid(n)

        ! Create field first
        obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec, &
                                     grid=LIS_obsvecGrid(n,k), &
                                     name="Observation"//vid(1)//vid(2), &
                                     rc=status)
        call LIS_verify(status, 'Field creation failed')

        ! After field creation
        write(LIS_logunit,*) '[DEBUG] Field created for nest:', n
        write(LIS_logunit,*) '[DEBUG] obs_ngrid value:', LIS_rc%obs_ngrid(n)

        ! Add to state first
        call ESMF_StateAdd(OBS_State(n), (/obsField(n)/), rc=status)
        call LIS_verify(status, 'State add failed')

        ! After state add
        write(LIS_logunit,*) '[DEBUG] Field added to state for nest:', n

        ! Now set attributes after field is created and added to state
        call ESMF_AttributeSet(obsField(n), "GRIDTYPE", "1D", rc=status)
        call LIS_verify(status, 'Error setting GRIDTYPE')

        ! Modified attribute setting
        call ESMF_AttributeSet(obsField(n), "Grid Number", &
             LIS_rc%obs_ngrid(n), itemCount=1, rc=status)
        call LIS_verify(status, 'Grid Number attribute failed for nest '//trim(adjustl(nest_str)))

        ! Set state attributes
        call ESMF_AttributeSet(OBS_State(n), "Number Of Observations", &
             LIS_rc%obs_ngrid(n), rc=status)
        call LIS_verify(status, 'Number Of Observations attribute failed')

        ! Read observation attributes
        write(LIS_logunit,*) '[INFO] Opening attributes for observations ', &
             trim(LIS_rc%obsattribfile(k))
        
        ! Add file existence check
        inquire(file=trim(LIS_rc%obsattribfile(k)), exist=file_exists)
        if (.not. file_exists) then
            write(LIS_logunit,*) '[ERR] Observation attributes file does not exist: ', &
                               trim(LIS_rc%obsattribfile(k))
            call LIS_endrun()
        endif
        
        ftn = LIS_getNextUnitNumber()
        open(ftn, file=trim(LIS_rc%obsattribfile(k)), status='old', iostat=status)
        
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Failed to open observation attributes file: ', &
                               trim(LIS_rc%obsattribfile(k))
            call LIS_endrun()
        endif
        
        ! Read the file with error checking
        read(ftn,*, iostat=status)
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Error reading line 1 from attributes file'
            call LIS_endrun()
        endif
        
        read(ftn,*, iostat=status) LIS_rc%nobtypes(k)
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Error reading nobtypes from attributes file'
            call LIS_endrun()
        endif
        
        write(LIS_logunit,*) '[DEBUG] Read nobtypes(k) =', LIS_rc%nobtypes(k)
        
        ! Safety check for nobtypes
        if (LIS_rc%nobtypes(k) <= 0 .or. LIS_rc%nobtypes(k) > 1000) then
            write(LIS_logunit,*) '[ERR] Invalid nobtypes value:', LIS_rc%nobtypes(k)
            call LIS_endrun()
        endif
        
        read(ftn,*, iostat=status)
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Error reading line 3 from attributes file'
            call LIS_endrun()
        endif

        ! Allocate arrays with error checking
        allocate(vname(LIS_rc%nobtypes(k)), stat=status)
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Failed to allocate vname array of size', LIS_rc%nobtypes(k)
            call LIS_endrun()
        endif
        
        allocate(varmax(LIS_rc%nobtypes(k)), stat=status)
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Failed to allocate varmax array of size', LIS_rc%nobtypes(k)
            call LIS_endrun()
        endif
        
        allocate(varmin(LIS_rc%nobtypes(k)), stat=status)
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Failed to allocate varmin array of size', LIS_rc%nobtypes(k)
            call LIS_endrun()
        endif

        do i=1, LIS_rc%nobtypes(k)
            read(ftn,fmt='(a40)') vname(i)
            read(ftn,*) varmin(i), varmax(i)
            write(LIS_logunit,*) '[INFO] ', vname(i), varmin(i), varmax(i)
        enddo
        call LIS_releaseUnitNumber(ftn)

        ! Allocate ssdev with error checking
        allocate(ssdev(LIS_rc%obs_ngrid(n)), stat=status)
        if (status /= 0) then
            write(LIS_logunit,*) '[ERR] Failed to allocate ssdev array of size', LIS_rc%obs_ngrid(n)
            call LIS_endrun()
        endif
        write(LIS_logunit,*) '[DEBUG] Allocated ssdev array with size =', size(ssdev)

        ! Handle perturbations
        if (trim(LIS_rc%perturb_obs(n)) /= "none") then
            write(LIS_logunit,*) '[DEBUG] Setting up perturbations for nest:', n
            
            ! Safe allocation with error checking
            allocate(obs_pert%vname(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%vname'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%perttype(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%perttype'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%ssdev(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%ssdev'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%stdmax(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%stdmax'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%zeromean(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%zeromean'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%tcorr(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%tcorr'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%xcorr(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%xcorr'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%ycorr(1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%ycorr'
                call LIS_endrun()
            endif
            
            allocate(obs_pert%ccorr(1,1), stat=status)
            if (status /= 0) then
                write(LIS_logunit,*) '[FATAL] Failed to allocate obs_pert%ccorr'
                call LIS_endrun()
            endif

            call LIS_readPertAttributes(1, LIS_rc%obspertAttribfile(k), obs_pert)
            ssdev = obs_pert%ssdev(1)

            pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec, grid=LIS_obsEnsOnGrid(n,k), &
                                           name="Observation"//vid(1)//vid(2), rc=status)
            call LIS_verify(status, 'Error ESMF_FieldCreate: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Created perturbation field for nest:', n

            ! Initialize the perturbations to zero
            call ESMF_FieldGet(pertField(n), localDE=0, farrayPtr=obs_temp, rc=status)
            call LIS_verify(status, 'Error ESMF_FieldGet: GPSdispObs_setup')
            obs_temp(:,:) = 0
            write(LIS_logunit,*) '[DEBUG] Initialized perturbations to zero'

            call ESMF_AttributeSet(pertField(n), "Perturbation Type", obs_pert%perttype(1), rc=status)
            call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Set Perturbation Type:', obs_pert%perttype(1)

            if (LIS_rc%obs_ngrid(n) > 0) then
                call ESMF_AttributeSet(pertField(n), "Standard Deviation", ssdev, &
                                      itemCount=LIS_rc%obs_ngrid(n), rc=status)
                call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
                write(LIS_logunit,*) '[DEBUG] Set Standard Deviation:', ssdev
            else
                write(LIS_logunit,*) '[WARN] No grid points for Standard Deviation'
            endif

            call ESMF_AttributeSet(pertField(n), "Std Normal Max", obs_pert%stdmax(1), rc=status)
            call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Set Std Normal Max:', obs_pert%stdmax(1)

            call ESMF_AttributeSet(pertField(n), "Ensure Zero Mean", obs_pert%zeromean(1), rc=status)
            call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Set Ensure Zero Mean:', obs_pert%zeromean(1)

            call ESMF_AttributeSet(pertField(n), "Temporal Correlation Scale", obs_pert%tcorr(1), rc=status)
            call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Set Temporal Correlation Scale:', obs_pert%tcorr(1)

            call ESMF_AttributeSet(pertField(n), "X Correlation Scale", obs_pert%xcorr(1), rc=status)
            call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Set X Correlation Scale:', obs_pert%xcorr(1)

            call ESMF_AttributeSet(pertField(n), "Y Correlation Scale", obs_pert%ycorr(1), rc=status)
            call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Set Y Correlation Scale:', obs_pert%ycorr(1)

            call ESMF_AttributeSet(pertField(n), "Cross Correlation Strength", obs_pert%ccorr(1,:), &
                                  itemCount=1, rc=status)
            call LIS_verify(status,'Error ESMF_AttributeSet: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Set Cross Correlation Strength'

            ! Add perturbation field to state
            call ESMF_StateAdd(OBS_Pert_State(n), (/pertField(n)/), rc=status)
            call LIS_verify(status,'Error ESMF_StateAdd: GPSdispObs_setup')
            write(LIS_logunit,*) '[DEBUG] Added perturbation field to OBS_Pert_State for nest:', n

            ! Deallocate perturbation structures
            deallocate(obs_pert%vname)
            deallocate(obs_pert%perttype)
            deallocate(obs_pert%ssdev)
            deallocate(obs_pert%stdmax)
            deallocate(obs_pert%zeromean)
            deallocate(obs_pert%tcorr)
            deallocate(obs_pert%xcorr)
            deallocate(obs_pert%ycorr)
            deallocate(obs_pert%ccorr)
            write(LIS_logunit,*) '[DEBUG] Deallocated perturbation structures for nest:', n

            ! Add missing grid number attribute to PERTURBATION field - only if pertField was created
            write(LIS_logunit,*) '[DEBUG] Setting Grid Number for nest',n,'dataset',k
            write(LIS_logunit,*) '[DEBUG] obs_ngrid size:',size(LIS_rc%obs_ngrid)
            write(LIS_logunit,*) '[DEBUG] obs_ngrid(',n,') =',LIS_rc%obs_ngrid(n)
            
            if (n > size(LIS_rc%obs_ngrid)) then
                write(LIS_logunit,*) '[ERR] Nest index out of bounds!',n,'>',size(LIS_rc%obs_ngrid)
                call LIS_endrun()
            endif
            
            if (LIS_rc%obs_ngrid(n) <= 0) then
                write(LIS_logunit,*) '[ERR] Invalid obs_ngrid for nest',n,'value:',LIS_rc%obs_ngrid(n)
                call LIS_endrun()
            endif
            
            write(ngrid_str,'(I0)') LIS_rc%obs_ngrid(n)
            call ESMF_AttributeSet(pertField(n),"Grid Number",&
                 LIS_rc%obs_ngrid(n),itemCount=LIS_rc%obs_ngrid(n),rc=status)
            call LIS_verify(status,'Error setting Grid Number in pertField: '//&
                             'nest='//trim(nest_str)//&
                             ' obs_ngrid='//trim(ngrid_str))

            ! Force attribute inheritance
            write(nest_str,'(I0)') n
            
            ! Get Grid Number from perturbation field
            call ESMF_AttributeGet(pertField(n), "Grid Number", gridNum, rc=status)
            call LIS_verify(status, 'Failed to get Grid Number from pertField for nest='//trim(nest_str))
            
            ! Set Grid Number on observation state
            call ESMF_AttributeSet(OBS_State(n), "Grid Number", gridNum, itemCount=1, rc=status)
            call LIS_verify(status,'Attribute copy failed nest='//trim(nest_str))
        endif

        ! Set GPS structure
        GPS_struc(n)%alarmhr = 3
        write(LIS_logunit,*) '[DEBUG] Set GPS alarm hour to:', GPS_struc(n)%alarmhr
        
        ! Cleanup
        deallocate(vname)
        deallocate(varmax)
        deallocate(varmin)
        deallocate(ssdev)
    enddo

  end subroutine GPSdispObs_setup

  subroutine gps_quality_control(Observations)
    type(GPS_observation), intent(inout) :: Observations(:)
    real, parameter :: GPS_DISP_THRESHOLD = 0.1 ! 10 cm
    integer :: i  ! Explicitly declare loop counter
    
    do i = 1, size(Observations)
      if (abs(Observations(i)%value) > GPS_DISP_THRESHOLD) then
        Observations(i)%assim = .false.
      endif
    enddo
  end subroutine

  ! Function to safely check if the Observations array is valid
  function is_observations_valid() result(valid)
    use LIS_logMod, only : LIS_logunit
    logical :: valid
    integer :: array_size
    integer :: status
    
    valid = .false.
    
    ! Basic check: is the array allocated?
    if (.not. allocated(Observations)) then
      write(LIS_logunit,*) '[WARN] Observations array is not allocated'
      return
    endif
    
    ! Get array size safely
    array_size = size(Observations)
    
    if (array_size <= 0) then
      write(LIS_logunit,*) '[WARN] Observations array has invalid size (≤ 0):', array_size
      return
    endif
    
    if (array_size > 1000000) then
      write(LIS_logunit,*) '[WARN] Observations array has suspiciously large size:', array_size
      return
    endif
    
    ! Array seems valid
    valid = .true.
    
  end function is_observations_valid

  ! Add a safe allocation function
  subroutine safely_allocate_observations(n)
    use LIS_logMod, only : LIS_logunit
    integer, intent(in) :: n
    integer :: status
    
    ! If already allocated, deallocate first
    if (allocated(Observations)) then
      deallocate(Observations, stat=status)
      if (status /= 0) then
        write(LIS_logunit,*) '[ERR] Failed to deallocate existing Observations array'
        call LIS_endrun()
      endif
    endif
    
    ! Sanity check on size
    if (n <= 0 .or. n > 1000000) then
      write(LIS_logunit,*) '[ERR] Invalid size for Observations array:', n
      call LIS_endrun()
    endif
    
    ! Allocate with error checking
    allocate(Observations(n), stat=status)
    if (status /= 0) then
      write(LIS_logunit,*) '[ERR] Failed to allocate Observations array of size', n
      call LIS_endrun()
    endif
    
    ! Initialize with safe defaults
    Observations%value = -9999.0
    Observations%error = -9999.0
    Observations%assim = .false.
    
  end subroutine safely_allocate_observations

  ! Add a safe deallocation function
  subroutine safely_deallocate_observations()
    use LIS_logMod, only : LIS_logunit
    integer :: status
    
    ! Check if already deallocated
    if (.not. allocated(Observations)) then
      write(LIS_logunit,*) '[INFO] Observations array already deallocated'
      return
    endif
    
    ! Deallocate with error checking
    deallocate(Observations, stat=status)
    if (status /= 0) then
      write(LIS_logunit,*) '[ERR] Failed to deallocate Observations array'
      call LIS_endrun()
    endif
    
  end subroutine safely_deallocate_observations

end module GPSdispObs_Mod
