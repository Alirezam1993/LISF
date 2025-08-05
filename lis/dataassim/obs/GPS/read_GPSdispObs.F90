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
! !ROUTINE: read_GPSdispObs
! \label{read_GPSdispObs}
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
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use GPSdispObs_Mod
  use LIS_fileIOMod
  use LIS_DAobservationsMod
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
!  Reads gridded GPS vertical displacement observations and packages 
!  them into an ESMF State with the necessary attributes for data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[k]    index of the observation type
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] perturbations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: dispfield 

  integer             :: iret
  real,    pointer    :: obsl(:)
  real                :: gpsgrid(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k))
  real                :: gpsgrid_glb(LIS_rc%obs_gnc(k), LIS_rc%obs_gnr(k))
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))

  character(len=LIS_CONST_PATH_LEN) :: GPSobsdir
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name
  logical             :: readflag
  logical             :: data_upd
  logical             :: data_upd_flag_local
  logical             :: data_upd_flag(LIS_npes)
  integer             :: fnd
  logical             :: alarmCheck
  integer             :: status, ftn, t, p, r, c, c1, r1, gid_map
  integer             :: valid_obs_count, valid_count_initial, valid_grid_count, valid_obs_set, invalid_index_count
  integer             :: obs_index  ! For direct array indexing fix
  real, allocatable   :: ssdev(:)   ! Standard deviation for DA system
  
  ! Initialize data update status
  data_upd = .false.
  
  ! Get data directory from state
  call ESMF_AttributeGet(OBS_State,"Data Directory", GPSobsdir, rc=status)
  call LIS_verify(status)
  
  ! Initially set data update status to false
  call ESMF_AttributeSet(OBS_State,"Data Update Status", .false., rc=status)
  call LIS_verify(status)
  
  ! Generate filename and check existence
  call GPS_filename(name, GPSobsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da)
  
  inquire(file=name, exist=file_exists)
  if(file_exists) then 
     call ESMF_AttributeSet(OBS_State,"File Status", .true., rc=status)
     call LIS_verify(status)
  else
     call ESMF_AttributeSet(OBS_State,"File Status", .false., rc=status)
     call LIS_verify(status)
  endif

  ! Check if we're at the correct time for assimilation
  ! GPS observations are available daily, unlike GRACE which is monthly
  if(LIS_rc%hr.eq.GPS_struc(n)%alarmHr .and. &
     LIS_rc%mn.eq.0 .and. &
     LIS_rc%ss.eq.0) then 
     alarmCheck = .true. 
  else
     alarmCheck = .false.
  endif

  ! Only proceed if we're at the right time and in the correct mode
  if(alarmCheck) then 
     if(LIS_rc%DAincrMode(n).eq.1) then 
        call GPS_filename(name, GPSobsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da)
        
        inquire(file=name, exist=file_exists)
        
        if(file_exists) then 
           readflag = .true. 
        else 
           readflag = .false.
        endif
        
        if (readflag) then 
           call ESMF_AttributeSet(OBS_State,"File Status", .true., rc=status)
           call LIS_verify(status)
           
           ! Set data averaging factor (1.0 for GPS, unlike monthly GRACE data)
           call ESMF_AttributeSet(OBS_State, name="Data averaging factor", &
                value=1.0, rc=status)
           call LIS_verify(status)

           write(LIS_logunit,*) '[INFO] Reading GPS data ', trim(name)

           ! Read GPS data - Fortran unformatted record with big-endian byte order
           ftn = LIS_getNextUnitNumber()
           open(ftn, file=trim(name), form='unformatted', status='old', &
                convert='big_endian', iostat=status)
           
           if(status /= 0) then
              write(LIS_logunit,*) '[WARN] Error opening GPS data file: ', trim(name)
              readflag = .false.
              call ESMF_AttributeSet(OBS_State, "File Status", .false., rc=status)
              call LIS_verify(status)
              call LIS_releaseUnitNumber(ftn)
              return
           endif
           
           ! Read GPS data as Fortran unformatted record
           read(ftn, iostat=status) gpsgrid_glb
           close(ftn)
           call LIS_releaseUnitNumber(ftn)
           
           if(status /= 0) then
              write(LIS_logunit,*) '[WARN] Error reading GPS data from file: ', trim(name)
              readflag = .false.
              call ESMF_AttributeSet(OBS_State, "File Status", .false., rc=status)
              call LIS_verify(status)
              return
           endif

           ! Debug: Check range of values read from file
           write(LIS_logunit,*) '[DEBUG] GPS data range - Min:', minval(gpsgrid_glb), ' Max:', maxval(gpsgrid_glb)
           write(LIS_logunit,*) '[DEBUG] GPS undefined count:', count(gpsgrid_glb == LIS_rc%udef)
           write(LIS_logunit,*) '[DEBUG] GPS -9999 count:', count(gpsgrid_glb == -9999.0)
           write(LIS_logunit,*) '[DEBUG] GPS global grid shape:', LIS_rc%obs_gnc(k), 'x', LIS_rc%obs_gnr(k)
           write(LIS_logunit,*) '[DEBUG] GPS -9999 count:', count(gpsgrid_glb == -9999.0)
           write(LIS_logunit,*) '[DEBUG] GPS global grid shape:', LIS_rc%obs_gnc(k), 'x', LIS_rc%obs_gnr(k)
           write(LIS_logunit,*) '[DEBUG] GPS local grid shape:', LIS_rc%obs_lnc(k), 'x', LIS_rc%obs_lnr(k)

           ! Extract local domain data
           gpsgrid = gpsgrid_glb( &
                LIS_ews_obs_halo_ind(n,LIS_localPet+1):&         
                LIS_ewe_obs_halo_ind(n,LIS_localPet+1), &
                LIS_nss_obs_halo_ind(n,LIS_localPet+1): &
                LIS_nse_obs_halo_ind(n,LIS_localPet+1))

           write(LIS_logunit,*) '[DEBUG] Local GPS data range - Min:', minval(gpsgrid), ' Max:', maxval(gpsgrid)
           write(LIS_logunit,*) '[DEBUG] Local GPS -9999 count:', count(gpsgrid == -9999.0)

           ! Initial check for any valid data in the file
           ! GPS data uses -9999.0 as undefined value
           fnd = 0 
           valid_count_initial = 0
           data_upd_flag_local = .false. 
           do r=1,LIS_rc%obs_lnr(k)
              do c=1,LIS_rc%obs_lnc(k)
                 ! Simplified QC: accept all data not equal to -9999.0 (like GRACE DA)
                 if(gpsgrid(c,r) /= -9999.0) then 
                    fnd = 1
                    valid_count_initial = valid_count_initial + 1
                 endif
              enddo
           enddo
           
           write(LIS_logunit,*) '[DEBUG] GPS initial QC found:', valid_count_initial, 'valid points, fnd=', fnd

           ! Get field and populate with observations
           call ESMF_StateGet(OBS_State, "Observation01", dispfield, rc=status)
           call LIS_verify(status, 'ESMF_StateGet failed in read_GPSdispObs')

           call ESMF_FieldGet(dispfield, localDE=0, farrayPtr=obsl, rc=status)
           call LIS_verify(status, 'ESMF_FieldGet failed in read_GPSdispObs')

           obsl(:) = -9999.0
           


           if(fnd.eq.0) then 
              obsl = LIS_rc%udef
           else
              valid_grid_count = 0
              valid_obs_set = 0
              invalid_index_count = 0
              

              ! FIXED: Use direct array indexing instead of broken gindex mapping
              obs_index = 0
              do r=1,LIS_rc%obs_lnr(k)
                 do c=1,LIS_rc%obs_lnc(k)
                    obs_index = obs_index + 1
                    ! Direct indexing approach - bypass broken gindex
                    if(obs_index <= LIS_rc%obs_ngrid(k)) then
                       valid_grid_count = valid_grid_count + 1
                       ! Simplified QC: like GRACE DA but accept negative values too
                       if(gpsgrid(c,r) /= -9999.0) then
                          obsl(obs_index) = gpsgrid(c,r)
                          valid_obs_set = valid_obs_set + 1
                       else
                          obsl(obs_index) = LIS_rc%udef
                       end if
                    else
                       invalid_index_count = invalid_index_count + 1
                    endif
                 enddo
              enddo
              
              write(LIS_logunit,*) '[DEBUG] Grid points with valid index:', valid_grid_count
              write(LIS_logunit,*) '[DEBUG] Grid points with invalid index (-1):', invalid_index_count
              write(LIS_logunit,*) '[DEBUG] Observations set in array:', valid_obs_set
           endif

           ! Recheck for valid observations after quality control
           fnd = 0
           valid_obs_count = 0
           ! FIXED: Use direct indexing for QC check as well
           obs_index = 0
           do r=1,LIS_rc%obs_lnr(k)
              do c=1,LIS_rc%obs_lnc(k)
                 obs_index = obs_index + 1
                 if(obs_index <= LIS_rc%obs_ngrid(k)) then
                    if(obsl(obs_index) /= LIS_rc%udef .and. &
                       obsl(obs_index) /= -9999.0) then 
                       fnd = 1
                       valid_obs_count = valid_obs_count + 1
                    endif
                 endif
              enddo
           enddo

           write(LIS_logunit,*) '[INFO] GPS data processing: Valid observations after QC = ', valid_obs_count

           ! Set data update flag
           if(fnd.eq.0) then 
              data_upd_flag_local = .false.
           else
              data_upd_flag_local = .true. 
           endif
           
           ! Gather update flags from all processes
#if(defined SPMD)
           call MPI_ALLGATHER(data_upd_flag_local, 1, &
                MPI_LOGICAL, data_upd_flag(:), &
                1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
           
           data_upd = .false.
           do p=1,LIS_npes
              data_upd = data_upd.or.data_upd_flag(p)
           enddo

           write(LIS_logunit,*) '[INFO] Read GPS data successfully'

           readflag = .false.

           if(data_upd) then              
              ! Set grid IDs and assimilation flags
              do t=1,LIS_rc%obs_ngrid(k)
                 gid(t) = t
                 if(obsl(t) /= -9999.0 .and. obsl(t) /= LIS_rc%udef) then 
                    assimflag(t) = 1
                 else
                    assimflag(t) = 0
                 endif
              enddo

              ! Mark data as updated
              call ESMF_AttributeSet(OBS_State, "Data Update Status", &
                   .true., rc=status)
              call LIS_verify(status)

              ! Set field attributes
              if(LIS_rc%obs_ngrid(k).gt.0) then 
                 call ESMF_AttributeSet(dispfield, "Grid Number", &
                      gid, itemCount=LIS_rc%obs_ngrid(k), rc=status)
                 call LIS_verify(status)

                 call ESMF_AttributeSet(dispfield, "Assimilation Flag", &
                      assimflag, itemCount=LIS_rc%obs_ngrid(k), rc=status)
                 call LIS_verify(status)
              endif
              
              ! Update Observations array for compatibility with other components
              valid_obs_count = 0
              do t=1,LIS_rc%obs_ngrid(k)
                 if (obsl(t) /= -9999.0 .and. obsl(t) /= LIS_rc%udef) then
                    Observations(t)%value = obsl(t)
                    Observations(t)%error = 20.0  ! GPS displacement error (from pertattribs file)
                    Observations(t)%assim = .true.
                    valid_obs_count = valid_obs_count + 1
                 else
                    Observations(t)%value = -9999.0
                    Observations(t)%error = -9999.0
                    Observations(t)%assim = .false.
                 endif
              enddo
              
              ! Set observation count attribute
              call ESMF_AttributeSet(OBS_State, "DA_Observation_Count", &
                   valid_obs_count, rc=status)
              call LIS_verify(status)
              
              ! GPS DA system uses Standard Deviation from pertattribs file
              ! This is set during module setup, not during observation reading
              ! The missing piece was that GPS never set this attribute in read function like GRACE does
              ! However, GPS should use the ssdev from setup perturbations, not set it here
              ! GRACE sets it here because it reads spatially variable errors from observation files
              ! GPS has uniform error from pertattribs file, already set in module setup
           endif
        else
           write(LIS_logunit,*) '[WARN] GPS data file not found: ', trim(name)
        endif
     endif
  endif
  
  
  
end subroutine read_GPSdispObs

!BOP
! 
! !ROUTINE: GPS_filename
! \label{GPS_filename}
! 
! !INTERFACE:  
subroutine GPS_filename(name, ndir, yr, mo, da)
! !USES:    
  implicit none
! !ARGUMENTS:     
  character(len=*)  :: name
  character(len=*)  :: ndir
  integer           :: yr, mo, da
!
! !DESCRIPTION:
!   
! This subroutine generates the standard filename for GPS displacement files.
!
!EOP
  character(len=4)  :: fyr
  character(len=2)  :: fmo, fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  
  ! Construct standard filename with simple path
  name = trim(ndir)//'/GPSdisp_'//trim(fyr)//trim(fmo)//trim(fda)//'0000.gs4r'
  
end subroutine GPS_filename
