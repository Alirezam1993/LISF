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
! !MODULE: noahmp401_gpsdisp_mod
! \label{noahmp401_gpsdisp_mod}
! 
! !DESCRIPTION:
! Module to hold shared variables for GPS displacement prediction
!
module noahmp401_gpsdisp_shared_mod
  implicit none
  private
  
  ! Module variables for NetCDF climatology file
  real, allocatable, save, public :: tws_climatology_grid(:)  ! Per observation point climatology
  logical, save, public          :: climatology_file_read = .false.
  real, save, public             :: tws_climatology = -9999.0      ! Store calculated climatology
  logical, save, public          :: climatology_calculated = .false.
  
end module noahmp401_gpsdisp_shared_mod

!BOP
! !ROUTINE: noahmp401_getgpsdisppred
! \label{noahmp401_getgpsdisppred}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
! Jan 2020: Jing Wang; Added GPS displacement functionality
! Oct 2024: Alireza Moghaddasi; Updated with Green's function
!
! !INTERFACE:
subroutine noahmp401_getgpsdisppred(n, k, obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use noahmp401_gpsdisp_DAlogMod
  use LIS_mpiMod
  use noahmp401_gpsdisp_shared_mod

!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))

!
! !DESCRIPTION:
!
!  Returns the GPS displacement obs pred (model's estimate of 
!  observations) for data assimilation.
!  Uses a Green's function approach to map water storage anomalies
!  to vertical displacement based on the Preliminary Reference Earth Model (PREM).
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[k] index of the observation type \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: t, i, m
  real                   :: tws(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real,dimension(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n)) :: tws_obs
  
    ! TWS climatological reference - UPDATED APPROACH
  real                   :: tws_mean_reference
  logical                :: use_dynamic_reference = .true.  ! Enable proper reference calculation
  
  logical                :: debug
  
  ! Variables for Green's function
  real, parameter        :: disk_radius_km = 14.0  ! Disk radius in km (for 25km grid)
  
  ! Forward declare external functions
  real, external         :: green_function_prem

  debug = .false.
  
      ! Read TWS climatological reference from NetCDF file
    if (use_dynamic_reference) then
        if (.not. climatology_calculated) then
            ! Read multi-year climatology from NetCDF file
            call read_tws_climatology_file(n, k)
            climatology_calculated = .true.
            write(LIS_logunit,*) '[GPS_DA] Successfully read TWS climatology from NetCDF file'
        endif
        ! tws_mean_reference will be applied per-pixel in observation space
    else
        ! Fallback to zero (original problematic approach)
        tws_mean_reference = 0.0
        write(LIS_logunit,*) '[GPS_DA] WARNING: Using zero TWS reference - may cause bias!'
    endif
  
  ! Calculate TWS directly from current model state
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
    tws(t) = &
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
    
    if (debug .and. t <= 3) then
      write(LIS_logunit,*) '[GPS_DA] Patch ', t, ': TWS=', tws(t), ' mm'
    endif
  enddo

  ! Convert patch space TWS to observation/ensemble space (absolute values)
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       tws,&
       tws_obs)
  
  ! Apply per-pixel climatology to create TWS anomalies in observation space
  if (use_dynamic_reference .and. allocated(tws_climatology_grid)) then
    do i=1,LIS_rc%obs_ngrid(k)
      do m=1,LIS_rc%nensem(n)
        if (tws_obs(i,m) /= -9999.0 .and. tws_climatology_grid(i) /= -9999.0) then
          ! Convert to anomaly using per-pixel climatology
          tws_obs(i,m) = tws_obs(i,m) - tws_climatology_grid(i)
          if (debug .and. i <= 3 .and. m == 1) then
            write(LIS_logunit,*) '[GPS_DA] Obs ', i, ': TWS_absolute=', tws_obs(i,m) + tws_climatology_grid(i), &
                                 ' Climatology=', tws_climatology_grid(i), ' Anomaly=', tws_obs(i,m)
          endif
        else
          if (tws_obs(i,m) /= -9999.0) then
            ! Fall back to domain average climatology
            tws_obs(i,m) = tws_obs(i,m) - tws_mean_reference
          endif
        endif
      enddo
    enddo
  else
    ! Fallback: use single reference value
    do i=1,LIS_rc%obs_ngrid(k)
      do m=1,LIS_rc%nensem(n)
        if (tws_obs(i,m) /= -9999.0) then
          tws_obs(i,m) = tws_obs(i,m) - tws_mean_reference
        endif
      enddo
    enddo
  endif

  ! Apply Green's function to map TWS anomaly to GPS displacement
  do i=1,LIS_rc%obs_ngrid(k)
     do m=1,LIS_rc%nensem(n)
      if (tws_obs(i,m) /= -9999.0) then
        ! Use Green's function to convert water mass to GPS displacement
        obs_pred(i,m) = green_function_prem(tws_obs(i,m), disk_radius_km)

        if (debug) then
          write(LIS_logunit,*) 'GPS point:', i, 'ensemble:', m
          write(LIS_logunit,*) 'TWS anomaly:', tws_obs(i,m), 'GPS displacement:', obs_pred(i,m)
        endif
      else
        obs_pred(i,m) = -9999.0
      endif
     enddo
  enddo

  if (debug) then
    write(LIS_logunit,*) 'GPS displacement prediction completed'
  endif

end subroutine noahmp401_getgpsdisppred

!BOP
! 
! !ROUTINE: green_function_prem
!  \label{green_function_prem}
! 
! !INTERFACE:   
real function green_function_prem(tws_anomaly, disk_radius_km) result(displacement)
! 
! !DESCRIPTION:    
!
! Green's function implementation based on the Preliminary Reference Earth Model (PREM)
! as described in Farrell (1972) and Wahr et al. (2013).
!
! This implementation calculates vertical displacement caused by water loading,
! using the formula:
!
! dr = Σ(l=0 to L) hl'*(4πGR/g(2l+1))*Pl(cos λ)
!
! where:
! - hl' is the elastic load Love number
! - G is Newton's gravitational constant
! - R is Earth's radius
! - g is gravitational acceleration at Earth's surface
! - Pl is the Legendre polynomial of degree l
! - λ is the angular distance from observation point to center of disk loading
!
! For a disk load, the formula is simplified by assuming the observation point
! is at the center of the disk (λ = 0), and we compute the vertical displacement
! directly from the TWS anomaly.
!
! inputs:
!  tws_anomaly - water mass anomaly in kg/m^2 (equivalent to mm of water)
!  disk_radius_km - radius of water loading disk in km
!
! outputs:
!  displacement - vertical displacement in mm
!
!EOP
    
  implicit none
    
  ! Input arguments
  real, intent(in) :: tws_anomaly         ! Water mass anomaly (kg/m^2)
  real, intent(in) :: disk_radius_km      ! Radius of water loading disk (km)
    
  ! Local variables
  real, parameter :: PI = 3.14159265359
  real :: earthRadius_km = 6371.0         ! Earth radius in km
  real :: disk_angular_radius             ! Angular radius of disk in radians
  real :: greens_coefficient              ! Green's function coefficient (mm per kg/m^2)
  
  ! Constants for disk load calculation
  real, parameter :: G = 6.67430e-11      ! Gravitational constant (m^3/kg/s^2)
  real, parameter :: gravity = 9.80665    ! Surface gravity (m/s^2)
  
  ! Elastic load Love numbers (from PREM model)
  ! These values are typically provided as tables. Here we use a single effective value
  ! that represents the weighted sum for a disk load of specified radius
  real, parameter :: effective_love_number = 0.6    ! Effective love number
  
  ! Calculate the angular radius of the disk
  disk_angular_radius = disk_radius_km / earthRadius_km
  
  ! The Green's function coefficient represents the vertical displacement
  ! per unit water mass loading for a disk of specified radius
  ! This is a simplification of the full sum in equation (4) in the literature
  
  ! For a disk load at the observation point (λ = 0),
  ! the coefficient approximately follows:
  greens_coefficient = (4.0 * PI * G * earthRadius_km * 1000.0 * effective_love_number) / &
                         (gravity * (1.0 - cos(disk_angular_radius)))
  
  ! Calculate displacement in mm
  ! Converting tws_anomaly from kg/m^2 to vertical displacement
  displacement = greens_coefficient * tws_anomaly
  
  ! Apply reasonableness limits to prevent extreme values
  if (abs(displacement) > 100.0) then
    displacement = sign(100.0, displacement)
  endif
    
  return
    
end function green_function_prem

!BOP
! 
! !ROUTINE: calculate_tws_climatology
!  \label{calculate_tws_climatology}
! 
! !INTERFACE:   
subroutine calculate_tws_climatology(n, tws_climatology)
! 
! !DESCRIPTION:    
!
! Calculates an approximate TWS climatology for use as reference in anomaly calculation.
! This is a temporary implementation that uses domain-average TWS from current state.
! 
! TODO: Replace with proper multi-year climatology calculation
!
! inputs:
!  n - nest index
!
! outputs:
!  tws_climatology - calculated climatological TWS value (mm)
!
!EOP
    
  use LIS_coreMod
  use LIS_logMod
  use LIS_constantsMod
  use noahmp401_lsmMod
    
  implicit none
    
  ! Input arguments
  integer, intent(in) :: n
  real, intent(out)   :: tws_climatology
    
  ! Local variables
  integer :: t, valid_patches
  real    :: tws_sum, tws_current
  
  tws_sum = 0.0
  valid_patches = 0
  
  ! Calculate domain-average TWS as approximation of climatology
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
    tws_current = &
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
    
    if (tws_current > 0.0) then  ! Basic validity check
      tws_sum = tws_sum + tws_current
      valid_patches = valid_patches + 1
    endif
  enddo
  
  if (valid_patches > 0) then
    tws_climatology = tws_sum / real(valid_patches)
  else
    ! Fallback value if no valid patches found
    tws_climatology = 1500.0  ! Reasonable estimate for many regions (mm)
    write(LIS_logunit,*) '[GPS_DA] WARNING: No valid TWS patches found, using fallback climatology'
  endif
  
  write(LIS_logunit,*) '[GPS_DA] Climatology calculated from ', valid_patches, ' patches'
  write(LIS_logunit,*) '[GPS_DA] Mean TWS: ', tws_climatology, ' mm'
  
  ! Sanity check - typical TWS ranges from 500-3000 mm
  if (tws_climatology < 100.0 .or. tws_climatology > 5000.0) then
    write(LIS_logunit,*) '[GPS_DA] WARNING: Unusual TWS climatology value: ', tws_climatology
  endif
    
end subroutine calculate_tws_climatology

!BOP
!
! !ROUTINE: read_tws_climatology_file
! \label{read_tws_climatology_file}
!
! !INTERFACE:
subroutine read_tws_climatology_file(n, k)
! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use netcdf
  use noahmp401_gpsdisp_shared_mod
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n  ! nest index
  integer, intent(in) :: k  ! observation type index
!
! !DESCRIPTION:
! This subroutine reads the TWS climatology from a NetCDF file
! and maps it to the observation grid points.
!
!EOP
  ! Local variables
  character(len=255) :: climatology_file
  integer :: ncid, varid, ndomain, ndomain_file
  real, allocatable :: lat1D(:), lon1D(:), tws_mean(:)
  integer :: i, ii
  real :: obs_lat, obs_lon, delta_d
  logical :: file_exist
  
  ! File path - read from LIS configuration 
  ! TODO: This should be read from LIS_rc configuration
  ! For now, use default relative path
  climatology_file = './TWS_climatology_CO_0p25_1D.nc'
  delta_d = 0.25  ! degree tolerance for spatial matching
  
  write(LIS_logunit,*) '[GPS_DA] Reading TWS climatology file: ', trim(climatology_file)
  
  ! Check if file exists
  inquire(file=climatology_file, exist=file_exist)
  if (.not. file_exist) then
    write(LIS_logunit,*) '[GPS_DA] ERROR: Climatology file does not exist: ', trim(climatology_file)
    call LIS_endrun()
  endif
  
  ! Open NetCDF file
  call LIS_verify(nf90_open(trim(climatology_file), NF90_NOWRITE, ncid), &
                  'Error opening TWS climatology NetCDF file')
  
  ! Read ndomain
  call LIS_verify(nf90_inq_varid(ncid, 'ndomain', varid), &
                  'nf90_inq_varid for ndomain failed')
  call LIS_verify(nf90_get_var(ncid, varid, ndomain_file), &
                  'nf90_get_var for ndomain failed')
  
  write(LIS_logunit,*) '[GPS_DA] Climatology file ndomain = ', ndomain_file
  
  ! Allocate arrays
  allocate(lat1D(ndomain_file))
  allocate(lon1D(ndomain_file))
  allocate(tws_mean(ndomain_file))
  
  ! Read coordinate arrays
  call LIS_verify(nf90_inq_varid(ncid, 'lat1D', varid), &
                  'nf90_inq_varid for lat1D failed')
  call LIS_verify(nf90_get_var(ncid, varid, lat1D), &
                  'nf90_get_var for lat1D failed')
  
  call LIS_verify(nf90_inq_varid(ncid, 'lon1D', varid), &
                  'nf90_inq_varid for lon1D failed')
  call LIS_verify(nf90_get_var(ncid, varid, lon1D), &
                  'nf90_get_var for lon1D failed')
  
  ! Read TWS mean array
  call LIS_verify(nf90_inq_varid(ncid, 'TWS_mean', varid), &
                  'nf90_inq_varid for TWS_mean failed')
  call LIS_verify(nf90_get_var(ncid, varid, tws_mean), &
                  'nf90_get_var for TWS_mean failed')
  
  ! Close NetCDF file
  call LIS_verify(nf90_close(ncid), 'Error closing TWS climatology file')
  
  ! Allocate observation grid climatology array
  if (.not. allocated(tws_climatology_grid)) then
    allocate(tws_climatology_grid(LIS_rc%obs_ngrid(k)))
  endif
  
  ! Map climatology to observation grid using nearest neighbor
  write(LIS_logunit,*) '[GPS_DA] Mapping climatology to ', LIS_rc%obs_ngrid(k), ' observation points'
  
  do i = 1, LIS_rc%obs_ngrid(k)
    obs_lat = LIS_obs_domain(n,k)%lat(i)
    obs_lon = LIS_obs_domain(n,k)%lon(i)
    
    ! Initialize with missing value
    tws_climatology_grid(i) = -9999.0
    
    ! Find nearest climatology point
    do ii = 1, ndomain_file
      if (lat1D(ii) /= -9999.0 .and. lon1D(ii) /= -9999.0 .and. &
          tws_mean(ii) /= -9999.0) then
        if (abs(lat1D(ii) - obs_lat) < delta_d .and. &
            abs(lon1D(ii) - obs_lon) < delta_d) then
          tws_climatology_grid(i) = tws_mean(ii)
          exit  ! Found match, move to next observation point
        endif
      endif
    enddo
  enddo
  
  ! Count valid mappings
  ndomain = count(tws_climatology_grid /= -9999.0)
  write(LIS_logunit,*) '[GPS_DA] Successfully mapped climatology to ', ndomain, &
                       ' out of ', LIS_rc%obs_ngrid(k), ' observation points'
  
  if (ndomain > 0) then
    ! Calculate statistics of mapped values
    tws_climatology = sum(tws_climatology_grid, mask=(tws_climatology_grid /= -9999.0)) / real(ndomain)
    write(LIS_logunit,*) '[GPS_DA] Climatology statistics: mean = ', tws_climatology, ' mm'
  else
    write(LIS_logunit,*) '[GPS_DA] WARNING: No valid climatology mappings found!'
  endif
  
  ! Clean up
  deallocate(lat1D, lon1D, tws_mean)
  
  write(LIS_logunit,*) '[GPS_DA] TWS climatology file reading completed'
  
end subroutine read_tws_climatology_file
