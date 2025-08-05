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
! !ROUTINE: write_GPSdispobs
! \label{write_GPSdispobs}
! 
! !REVISION HISTORY: 
! 25Jan2008: Sujay Kumar; Initial Specification
! 28Sep2011: Ben Zaitchik; modified for GPSdisp
!
! !INTERFACE: 
subroutine write_GPSdispobs(n, k, OBS_State)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_historyMod
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none

! !ARGUMENTS: 

  integer,     intent(in)  :: n 
  integer,     intent(in)  :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION: 
! 
! writes the transformed (interpolated/upscaled/reprojected)  
! GPS displacement observations to a file
! 
!EOP
  type(ESMF_Field)         :: dispfield
  logical                  :: data_update
  real, pointer            :: dispobs(:)
  real                     :: dispobs_unsc(LIS_rc%obs_ngrid(k))
  character(len=LIS_CONST_PATH_LEN) :: obsname
  integer                  :: ftn
  integer                  :: status

  call ESMF_AttributeGet(OBS_State, "Data Update Status", & 
       data_update, rc=status)
  call LIS_verify(status)

  if(data_update) then 
     
     call ESMF_StateGet(OBS_State, "Observation01", dispfield, &
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(dispfield, localDE=0, farrayPtr=dispobs, rc=status)
     call LIS_verify(status)

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call GPS_dispobsname(obsname)        

        call LIS_create_output_directory('DAOBS')
        open(ftn, file=trim(obsname), form='unformatted')
     endif

     call LIS_writevar_gridded_obs(ftn, n, k, dispobs)
     
     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif

  endif  

end subroutine write_GPSdispobs

!BOP
! !ROUTINE: GPS_dispobsname
! \label{GPS_dispobsname}
! 
! !INTERFACE: 
subroutine GPS_dispobsname(obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)      :: obsname
! 
! !DESCRIPTION: 
! This routine generates the filename for GPS displacement observations output.
! 
!EOP

  character(len=6) :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2)') &
       LIS_rc%yr, LIS_rc%mo

  obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//'/'//cdate1//   &
            '.1gs4r'
  
end subroutine GPS_dispobsname
