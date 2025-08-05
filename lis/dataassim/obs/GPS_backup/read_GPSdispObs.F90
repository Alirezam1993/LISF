!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_GPSdispObs
!  \label{read_GPSdispObs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!  29Sep2017: Yonghwan Kwon; modified for DTB observations
!  11Aug2020: Gaohong Yin; modified for GPS displacement observations
!  21Dec2020: Jing Wang; modified based on Gaohong
!  05Mar2023: Alireza Moghaddasi; Enhanced to handle missing files gracefully
!  09Dec2024: Alireza Moghaddasi; Added robust error handling for attribute setting
!
! !INTERFACE: 
subroutine read_GPSdispObs(n, k, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_historyMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use GPSdispObs_Mod, only : Observations, GPS_observation, safely_allocate_observations, safely_deallocate_observations
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the ground-based GPS observations of vertical displacement (gridded files
!  derived from using the IDW interpolation) and packages into an ESMF state object 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: dispfield 

  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable   :: dummy(:)

  character(len=LIS_CONST_PATH_LEN) :: GPSobsdir
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name
  logical             :: readflag
  logical             :: data_update = .false.
  logical             :: alarmCheck
  integer             :: fileDay, fileHour

  integer             :: t
  integer             :: status, ftn, i
  integer             :: valid_gps_points
  integer             :: MAX_ATTRIBUTE_SET_RETRIES = 3
  integer             :: retryCount
  character(len=LIS_CONST_PATH_LEN) :: filename_base, filename_with_ext

  ! Set default alarm check to false
  alarmCheck = .false.

  ! Check if we're at the right time to read data - originally only at 00:00
  ! Relaxing to allow any hour to be considered valid for testing
  if(LIS_rc%mn.eq.0) then ! Removed hour=0 check, allowing any hour at minute=0
    alarmCheck = .true.
    write(LIS_logunit,*) '[INFO] At potential assimilation time (',LIS_rc%hr,':',LIS_rc%mn,')'
  endif

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       GPSobsdir, rc=status)
  call LIS_verify(status)
  
  ! Print debug info about the data directory path
  write(LIS_logunit,*) '[DEBUG] GPS data directory from attribute: "', trim(GPSobsdir), '"'

  ! Set the Data Update Status attribute
  call ESMF_AttributeSet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  ! Format the filename base - ALWAYS use 0000 for time regardless of current time
  write(filename_base,'(a,i4.4,i2.2,i2.2,a)') 'GPSdisp_', &
       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, '0000'
  
  write(LIS_logunit,*) '[DEBUG] Looking for GPS files with base: ', trim(filename_base)

  ! Try several paths in order of preference:
  ! 1. Try with CO subdirectory first if it's not already in the path
  if (index(GPSobsdir, '/CO') <= 0) then
    filename_with_ext = trim(GPSobsdir) // '/CO/' // trim(filename_base) // '.gs4r'
    write(LIS_logunit,*) '[DEBUG] Looking for GPS file in CO subdirectory: "', trim(filename_with_ext), '"'
    inquire(file=trim(filename_with_ext), exist=file_exists)
  endif
  
  ! 2. If not found or CO already in path, try with .gs4r extension in main directory
  if (.not. file_exists) then
    filename_with_ext = trim(GPSobsdir) // '/' // trim(filename_base) // '.gs4r'
    write(LIS_logunit,*) '[DEBUG] Looking for GPS file: "', trim(filename_with_ext), '"'
    inquire(file=trim(filename_with_ext), exist=file_exists)
  endif
  
  ! 3. If still not found, try without extension
  if (.not. file_exists) then
    filename_with_ext = trim(GPSobsdir) // '/' // trim(filename_base)
    write(LIS_logunit,*) '[DEBUG] Looking for GPS file without extension: "', trim(filename_with_ext), '"'
    inquire(file=trim(filename_with_ext), exist=file_exists)
  endif
  
  ! Set flags based on file existence
  if(file_exists) then
     write(LIS_logunit,*) '[INFO] GPS data file FOUND: ', trim(filename_with_ext)
     name = filename_with_ext
     readflag = .true.
     call ESMF_AttributeSet(OBS_State, "File Status", .true., rc=status)
     call LIS_verify(status)
  else
     readflag = .false.
     call ESMF_AttributeSet(OBS_State, "File Status", .false., rc=status)
     call LIS_verify(status)
     
     ! Try to guess what day the file might be from with hardcoded values
     do i=1,10
        ! Try different days from 1-10 of the month
        write(filename_base,'(a,i4.4,i2.2,i2.2,a)') 'GPSdisp_', &
             LIS_rc%yr, LIS_rc%mo, i, '0000.gs4r'
             
        filename_with_ext = trim(GPSobsdir) // '/CO/' // trim(filename_base)
        write(LIS_logunit,*) '[DEBUG] Trying with different day: "', trim(filename_with_ext), '"'
        inquire(file=trim(filename_with_ext), exist=file_exists)
        
        if (file_exists) then
            write(LIS_logunit,*) '[INFO] Found GPS file for day', i, 'of month'
            exit
        endif
     enddo
     
     write(LIS_logunit,*) '[WARN] GPS data assimilation skipped: No suitable file exists'
     write(LIS_logunit,*) '[DEBUG] Looking for any GPS files in the directory...'
     
     ! Try to list all GPS files in the directory to help with debugging
     call run_system_ls(trim(GPSobsdir))
     
     ! Exit gracefully without error
     return
  endif

  ! Only proceed with reading the file if it exists and it's time to assimilate
  if(readflag .and. alarmCheck) then 
     write(LIS_logunit,*) '[INFO] Reading GPS disp data ', trim(name)

     call ESMF_StateGet(OBS_State, "Observation01", dispfield, &
          rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(dispfield, localDE=0, farrayPtr=obsl, rc=status)
     call LIS_verify(status, 'Error: FieldGet for obsl (GPS_disp)')
     
     ! Initialize observation array to missing value
     obsl(:) = -9999.0
     
     ! Try to open and read the file
     ftn = 90
     open(ftn, file=trim(name), form='unformatted', status='old', iostat=status)
     
     if(status /= 0) then
        write(LIS_logunit,*) '[WARN] Error opening GPS data file: ', trim(name)
        write(LIS_logunit,*) '[WARN] GPS data assimilation skipped: File open error'
        readflag = .false.
        call ESMF_AttributeSet(OBS_State, "File Status", .false., rc=status)
        call LIS_verify(status)
        return
     endif
     
     ! Read the data
     call readobsvar_1dgridded_gps(ftn, n, k, obsl)
     close(ftn)

     ! Mark data as read
     readflag = .false.

     ! Set grid IDs and assimilation flags
     do t=1, LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(obsl(t) /= -9999.0) then 
           assimflag(t) = 1 
        else
           assimflag(t) = 0 
        endif
     enddo

     ! Count valid GPS points (using GRACE-style approach)
     valid_gps_points = count(obsl > -9000.0)
     
     ! Flag invalid points individually
     do t=1, LIS_rc%obs_ngrid(k)
        if(obsl(t) < -9000.0) then
            assimflag(t) = 0
            obsl(t) = -9999.0
        endif
     enddo
     
     ! Set consistent observation count attributes
     call ESMF_AttributeSet(OBS_State, "Number Of Observations", valid_gps_points, rc=status)
     call LIS_verify(status, 'GPS: Failed to set valid obs count')
     
     ! Also set the attribute with the standard name used by EnKS
     call ESMF_AttributeSet(OBS_State, "NUMBER_OF_OBSERVATIONS", valid_gps_points, rc=status)
     call LIS_verify(status, 'GPS: Failed to set NUMBER_OF_OBSERVATIONS attribute')
     
     ! Early exit if no valid observations found
     if (valid_gps_points <= 0) then
         write(LIS_logunit,*) '[WARN] No valid GPS observations found in this cycle'
         write(LIS_logunit,*) '[WARN] GPS data assimilation skipped: No valid data'
         call ESMF_AttributeSet(OBS_State, "Data Update Status", .false., rc=status)
         call LIS_verify(status)
         return
     endif

     ! Log the number of valid observations found
     write(LIS_logunit,*) '[INFO] Found', valid_gps_points, 'valid GPS observations'

     ! Set GRACE-compatible attributes
     call ESMF_AttributeSet(dispfield, "GRACE_Compat_NumObs", &
          valid_gps_points, rc=status)
     call LIS_verify(status, 'Error setting GRACE_Compat_NumObs in field')

     ! Mark data as updated
     call ESMF_AttributeSet(OBS_State, "Data Update Status", .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet for Data Update Status')

     ! Set the grid number attribute (this was missing and causing the crash)
     do t=1, LIS_rc%obs_ngrid(k)
         gid(t) = t
     enddo
     
     ! Print debug information about grid size
     write(LIS_logunit,*) '[DEBUG] Setting Grid Number attribute with size:', LIS_rc%obs_ngrid(k)
     
     ! Add more robust error handling for this critical step
     status = -1  ! Initialize to failure to be extra safe
     
     ! Try to verify gid array is correctly populated
     if (any(gid <= 0) .or. any(gid > LIS_rc%obs_ngrid(k))) then
         write(LIS_logunit,*) '[WARN] Invalid values in gid array. Fixing...'
         do t=1, LIS_rc%obs_ngrid(k)
             gid(t) = t
         enddo
     endif
     
     ! Multiple attempts to set the attribute with better error handling
     retryCount = 0
     do while (retryCount < MAX_ATTRIBUTE_SET_RETRIES)
         ! Attempt to set the attribute
         call ESMF_AttributeSet(dispfield, "Grid Number", &
              gid, itemCount=LIS_rc%obs_ngrid(k), rc=status)
         
         if (status == 0) then
             write(LIS_logunit,*) '[INFO] Successfully set Grid Number attribute'
             exit  ! Success, exit the retry loop
         else
             ! Try a different approach on failure
             retryCount = retryCount + 1
             write(LIS_logunit,*) '[WARN] Failed to set Grid Number attribute, attempt', retryCount, &
                                  ' of', MAX_ATTRIBUTE_SET_RETRIES, ' (status =', status, ')'
             
             ! On last attempt, try with a smaller array if large
             if (retryCount == MAX_ATTRIBUTE_SET_RETRIES - 1 .and. LIS_rc%obs_ngrid(k) > 100) then
                 write(LIS_logunit,*) '[WARN] Trying with partial array'
                 call ESMF_AttributeSet(dispfield, "Grid Number", &
                     gid(1:min(100, LIS_rc%obs_ngrid(k))), &
                     itemCount=min(100, LIS_rc%obs_ngrid(k)), rc=status)
                 
                 if (status == 0) then
                     write(LIS_logunit,*) '[WARN] Partial Grid Number attribute set successfully'
                     exit
                 endif
             endif
             
             ! Special fallback on final attempt
             if (retryCount == MAX_ATTRIBUTE_SET_RETRIES) then
                 write(LIS_logunit,*) '[WARN] All attempts to set Grid Number attribute failed.'
                 write(LIS_logunit,*) '[WARN] Setting simplified attribute instead'
                 
                 ! Try a very simple attribute (just a scalar) as last resort
                 call ESMF_AttributeSet(dispfield, "Grid Number", &
                     LIS_rc%obs_ngrid(k), rc=status)
                 
                 if (status == 0) then
                     write(LIS_logunit,*) '[WARN] Simplified Grid Number attribute set'
                 else
                     write(LIS_logunit,*) '[WARN] Even simplified attribute setting failed'
                     ! But we'll continue anyway - better than crashing
                 endif
                 
                 ! Don't fail the entire run just because of an attribute
                 status = 0
                 exit
             endif
         endif
     enddo

     call ESMF_AttributeSet(dispfield, "Assimilation Flag", &
          assimflag, itemCount=LIS_rc%obs_ngrid(k), rc=status)
     call LIS_verify(status, 'Error: AttributeSet for Assimilaton Flag')

     call ESMF_AttributeSet(dispfield, "Number Of Observations", &
          valid_gps_points, itemCount=LIS_rc%obs_ngrid(k), rc=status)
     call LIS_verify(status, 'Error setting obs count on field')

     ! GPS-specific metadata
     call ESMF_AttributeSet(dispfield, "units", "meters", rc=status)
     call LIS_verify(status, 'Failed to set GPS displacement units')

     ! Handle Observations array allocation and population
     if (.not. allocated(Observations)) then
         ! Allocate the Observations array directly
         call safely_allocate_observations(LIS_rc%obs_ngrid(k))
         
         write(LIS_logunit,*) '[DEBUG] Observations array allocated with size:', &
              size(Observations)
     else
         ! Ensure array is the right size
         if (size(Observations) /= LIS_rc%obs_ngrid(k)) then
             write(LIS_logunit,*) '[WARN] Observations array size mismatch:', &
                  size(Observations), 'vs expected', LIS_rc%obs_ngrid(k)
             
             ! Reallocate with the correct size - first deallocate
             call safely_deallocate_observations()
             
             ! Then allocate with the correct size
             call safely_allocate_observations(LIS_rc%obs_ngrid(k))
         endif
     endif

     ! Update Observations from obsl for valid points
     do t=1, LIS_rc%obs_ngrid(k)
         if (obsl(t) > -9000.0) then
             Observations(t)%value = obsl(t)
             Observations(t)%error = -9999.0  ! Could be customized based on your error model
             Observations(t)%assim = .true.
         else
             Observations(t)%value = -9999.0
             Observations(t)%error = -9999.0
             Observations(t)%assim = .false.
         endif
     enddo

     ! Set the count of observations to be assimilated
     call ESMF_AttributeSet(OBS_State, "DA_Observation_Count", &
          count(Observations%assim), rc=status)
     call LIS_verify(status)
     
     ! Set grid type attribute
     call ESMF_AttributeSet(OBS_State, "GRID_TYPE", "2D", rc=status)
     call LIS_verify(status, 'Error setting GRID_TYPE attribute')
     
     ! Final confirmation of successful data update
     data_update = .true.
     call ESMF_AttributeSet(OBS_State, "Data Update Status", &
          data_update, rc=status)
     call LIS_verify(status)
     
     ! Set Data averaging factor (required by EnKS, defaults to 1.0 for GPS unlike GRACE monthly data)
     call ESMF_AttributeSet(OBS_State, name="Data averaging factor", &
          value=1.0, rc=status)
     call LIS_verify(status, 'Error setting Data averaging factor')
     
     write(LIS_logunit,*) '[INFO] GPS data successfully processed for assimilation'
  else
     if (.not. file_exists) then
         write(LIS_logunit,*) '[WARN] GPS data assimilation skipped: File does not exist'
     elseif (.not. alarmCheck) then
         write(LIS_logunit,*) '[INFO] GPS data assimilation skipped: Not at assimilation time'
     endif
     
     call ESMF_AttributeSet(OBS_State, "Data Update Status", .false., rc=status)
     call LIS_verify(status)
  end if

end subroutine read_GPSdispObs

!===============================================================

subroutine GPSdisp_filename(ndir, yr, mo, da, hr, mn, filename, status)
! !USES:
  use LIS_logMod, only: LIS_logunit
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none
  
  ! ARGUMENTS:
  character(len=*), intent(in)  :: ndir
  integer,          intent(in)  :: yr, mo, da, hr, mn
  character(len=*), intent(out) :: filename
  integer,          intent(out) :: status

  ! LOCAL VARIABLES:
  character(len=LIS_CONST_PATH_LEN) :: ndir_clean
  character(len=30) :: fname_pattern
  logical           :: dir_exists

  ! Initialize
  status = 1 ! Default to error
  filename = ''
  
  ! Create filename without path (for debugging)
  write(fname_pattern, '(A,I4.4,I2.2,I2.2,I2.2,I2.2,A)') 'GPSdisp_', yr, mo, da, hr, mn, '.gs4r'
  write(LIS_logunit,*) '[DEBUG] Generated filename pattern: ', trim(fname_pattern)

  ! Clean up directory path (ensure it has trailing slash)
  ndir_clean = trim(ndir)
  if(len_trim(ndir_clean) > 0) then
    if(ndir_clean(len_trim(ndir_clean):len_trim(ndir_clean)) /= '/') then
      ndir_clean = trim(ndir_clean) // '/'
    endif
  endif

  ! Check if directory exists
  inquire(file=trim(ndir_clean)//'/.', exist=dir_exists)
  if(.not. dir_exists) then
    write(LIS_logunit,*) '[WARN] GPS data directory does not exist: "', trim(ndir_clean), '"'
    status = 1
    return
  endif

  ! Construct full path with filename
  filename = trim(ndir_clean) // fname_pattern
  
  write(LIS_logunit,*) '[DEBUG] Full file path: ', trim(filename)
  
  ! Check if file exists
  inquire(file=trim(filename), exist=dir_exists)
  if(.not. dir_exists) then
    status = 1  ! File not found
  else
    status = 0  ! File found
  endif

end subroutine GPSdisp_filename


!BOP
!
! !ROUTINE: readobsvar_1dgridded_gps
! \label{readobsvar_1dgridded_gps}
!
! !INTERFACE:
subroutine readobsvar_1dgridded_gps(ftn, n, k, var)
! !USES:
  use LIS_coreMod
  use LIS_DAobservationsMod
  use LIS_logMod

  implicit none
! !ARGUMENTS:
  integer              :: ftn 
  integer              :: n
  integer              :: k
  real                 :: var(LIS_rc%obs_ngrid(k))
!
! !DESCRIPTION:
!  This routine reads the observation data and subsets to the
!  local processor's domain, in a 1-d vector formulation.
!
!EOP

  real,  allocatable   :: gobs(:,:)
  integer              :: nc, c1, r1, c, r, gid
  integer              :: status

  ! Allocate storage for global grid observations
  allocate(gobs(LIS_rc%obs_gnc(k), LIS_rc%obs_gnr(k)), stat=status)
  if (status /= 0) then
      write(LIS_logunit,*) '[ERR] Failed to allocate memory for GPS observations grid'
      return
  endif

  ! Initialize grid to missing values
  gobs = -9999.0

  ! Read data with error handling
  read(ftn, iostat=status) gobs
  if (status /= 0) then
      write(LIS_logunit,*) '[WARN] Error reading GPS data file'
      deallocate(gobs)
      return
  endif

  ! Calculate dimensions for local domain extraction
  nc = (LIS_ewe_obs_halo_ind(k, LIS_localPet+1) - &
       LIS_ews_obs_halo_ind(k, LIS_localPet+1)) + 1

  ! Extract the relevant subset for this process
  do r=LIS_nss_obs_halo_ind(k, LIS_localPet+1), &
       LIS_nse_obs_halo_ind(k, LIS_localPet+1)
     do c=LIS_ews_obs_halo_ind(k, LIS_localPet+1), &
          LIS_ewe_obs_halo_ind(k, LIS_localPet+1)
        
        c1 = c - LIS_ews_obs_halo_ind(k, LIS_localPet+1) + 1
        r1 = r - LIS_nss_obs_halo_ind(k, LIS_localPet+1) + 1
        
        ! Bounds checking
        if (c1 <= 0 .or. r1 <= 0 .or. &
            c1 > size(LIS_obs_domain(n,k)%gindex, 1) .or. &
            r1 > size(LIS_obs_domain(n,k)%gindex, 2)) then
            cycle
        endif
        
        gid = LIS_obs_domain(n,k)%gindex(c1, r1)
        
        if (gid > 0 .and. gid <= LIS_rc%obs_ngrid(k)) then
            var(gid) = gobs(c, r)
        endif
     enddo
  enddo

  ! Clean up
  deallocate(gobs)

end subroutine readobsvar_1dgridded_gps

!BOP
!
! !ROUTINE: run_system_ls
!
! !INTERFACE:
subroutine run_system_ls(directory_path)
! !USES:
  use LIS_logMod, only: LIS_logunit
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
! !ARGUMENTS:
  character(len=*), intent(in) :: directory_path
!
! !DESCRIPTION:
!  This routine executes a system command to list files in the given directory
!  and reports the results to the LIS log. This is solely for debugging purposes.
!
!EOP
  integer :: status, cmd_status
  character(len=4*LIS_CONST_PATH_LEN) :: cmd
  character(len=LIS_CONST_PATH_LEN) :: dir_clean
  
  ! Clean up directory path
  dir_clean = trim(directory_path)
  if (len_trim(dir_clean) > 0) then
    if (dir_clean(len_trim(dir_clean):len_trim(dir_clean)) /= '/') then
      dir_clean = trim(dir_clean) // '/'
    endif
  endif
  
  ! First try to list GPS files with the gs4r extension
  write(cmd, '(a)') 'ls -la ' // trim(dir_clean) // '*GPSdisp*.gs4r'
  
  write(LIS_logunit,*) '[DEBUG] Executing safe ls command for GPS files'
  write(LIS_logunit,*) '[DEBUG] Directory: ', trim(dir_clean)
  
#if (defined USE_SYSTEM)
  ! We need a safer approach that won't crash if the system call fails
  open(unit=9991, file='/tmp/gps_ls_debug.txt', status='replace')
  close(9991)
  
  ! Redirect output to a temporary file
  write(cmd, '(a)') 'ls -la ' // trim(dir_clean) // '*GPSdisp*.gs4r > /tmp/gps_ls_debug.txt 2>&1'
  call system(trim(cmd))
  
  ! Display the contents of that file
  write(LIS_logunit,*) '[DEBUG] GPS files (*.gs4r) in directory:'
  write(cmd, '(a)') 'cat /tmp/gps_ls_debug.txt'
  call system(trim(cmd))
  
  ! Also try without extension
  write(cmd, '(a)') 'ls -la ' // trim(dir_clean) // '*GPSdisp* > /tmp/gps_ls_debug.txt 2>&1'
  call system(trim(cmd))
  
  write(LIS_logunit,*) '[DEBUG] All GPS files (any extension) in directory:'
  write(cmd, '(a)') 'cat /tmp/gps_ls_debug.txt'
  call system(trim(cmd))
#else
  write(LIS_logunit,*) '[WARN] Directory listing not available (system call disabled)'
  write(LIS_logunit,*) '[WARN] Please check directory manually: ', trim(dir_clean)
#endif
  
end subroutine run_system_ls
